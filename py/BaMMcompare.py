'''
This script is used for comparing a given motif to motifs in a given database.
Prerequisite: all motifs must be in BaMM format.
'''

import argparse
from multiprocessing import Pool
import numpy as np
import logging
import sys

from utils import calculate_H_model_bg, calculate_H_model, model_sim, update_models, filter_pwms, parse_meme, write_meme


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="input MEME-format file with multiple PWMs")
    parser.add_argument('output_file', help="name the output file with path")

    parser.add_argument('--model_db', default=None, help="specify the path to database that you want to search for")

    parser.add_argument('--db_order', type=int, default=4, help="the order of motifs in the database. Default: 4")
    parser.add_argument('--query_order', type=int, default=0, help="the order of query motif. Default: 0")

    parser.add_argument('--n_neg_perm', type=int, default=10)
    parser.add_argument('--highscore_fraction', type=float, default=0.1)
    parser.add_argument('--evalue_threshold', type=float, default=0.1)
    parser.add_argument('--pvalue_threshold', type=float, default=0.01)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--min_overlap', type=int, default=2, help="minimal overlaps between PWMs. Default: 2")
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('--output_score_file', default=None)

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    query_file = args.input_file
    db_file = args.model_db
    out_file = args.output_file
    output_score_file = args.output_score_file

    min_overlap = args.min_overlap


    # print logs out
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger_fmt = '%(asctime)s [%(levelname)s]  %(message)s'
    formatter = logging.Formatter(logger_fmt)
    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logger.setLevel(logging.INFO)

    # parse input query meme file
    model_set = parse_meme(query_file)
    
    # pre-compute entropy for all the models
    models = update_models(model_set['models'])

    db_models = models
    # if model_db is not given, search model files against themselves
    if db_file != None:
        model_db = parse_meme(db_file)
        db_models = update_models(model_db['models'])
    db_size = len(db_models)

    #filter models using affinity propagation
    new_models = filter_pwms(models, min_overlap)
    # update the model set after filtering
    new_model_set = {}
    new_model_set['version']    = model_set['version']
    new_model_set['alphabet']   = model_set['alphabet']
    new_model_set['bg_freq']    = model_set['bg_freq']
    new_model_set['models']     = new_models
    # write out the meme file
    write_meme(new_model_set, out_file)

    if output_score_file != None:
        def init_workers():
            global highscore_fraction_g
            highscore_fraction_g = args.highscore_fraction
            global evalue_thresh_g
            evalue_thresh_g = args.evalue_threshold
            np.random.seed(args.seed)
            global db_models_g
            db_models_g = db_models
            global db_size_g
            db_size_g = db_size
            global n_neg_perm_g
            n_neg_perm_g = args.n_neg_perm
            global min_overlap_g
            min_overlap_g = args.min_overlap

        logger.info('Queuing %s search jobs', len(models))

        with open(output_score_file, 'w') as out:
            print('model_id', 'db_id', 'simscore', 'e-value',
                  'start_query', 'end_query', 'start_hit', 'end_hit', 'bg_score', 'cross_score', 'pad_score',
                  sep='\t', file=out)
            with Pool(args.n_processes, initializer=init_workers) as pool:
                jobs = []
                for model in models:
                    job = pool.apply_async(motif_search, args=(model,))
                    jobs.append(job)

                total_jobs = len(jobs)
                for job_index, job in enumerate(jobs, start=1):
                    hits = job.get()
                    hits.sort(key=lambda x: x[3])
                    for hit in hits:
                        print(*hit, sep='\t', file=out)
                        logger.info('Finished (%s/%s)', job_index, total_jobs)


def motif_search(model):
    model_pwm   = model['pwm']
    model_len   = len(model_pwm)
    bg_freq     = model['bg_freq']

    assert model_len > 1

    # use shuffled pwms to estimate the p-value under the null.
    shuffled_dists = []
    for _ in range(n_neg_perm_g):
        # create shuffled model
        shuffle_model = {}
        # calculate a locality preserving permutation
        Z = np.random.normal(0, 1, model_len)
        shuffle_ind = np.argsort(2 * Z - np.arange(model_len))
        shuffle_model['pwm'] = model_pwm[shuffle_ind]

        # nucleotide substitution
        for j in range(model_len):
            # 1) switch A->T with p=0.5 per column
            if np.random.normal(0, 1) > 0.5:
                shuffle_model['pwm'][j][0], shuffle_model['pwm'][j][3] = \
                    shuffle_model['pwm'][j][3], shuffle_model['pwm'][j][0]

            # 2) switch C->G with p=0.5 per column
            if np.random.normal(0, 1) > 0.5:
                shuffle_model['pwm'][j][1], shuffle_model['pwm'][j][2] = \
                    shuffle_model['pwm'][j][2], shuffle_model['pwm'][j][1]

        shuffle_model['H_model']    = calculate_H_model(shuffle_model['pwm'])
        shuffle_model['H_model_bg'] = calculate_H_model_bg(shuffle_model['pwm'], bg_freq)

        for db_model in db_models_g:
            suffle_sim, *_ = model_sim(shuffle_model, db_model, min_overlap_g)
            shuffled_dists.append(suffle_sim)

    # we are fitting only the tail of the null scores with an exponential
    # distribution
    sorted_null = np.sort(shuffled_dists)
    N_neg = len(sorted_null)
    high_scores = sorted_null[-int(N_neg * highscore_fraction_g):]
    high_score = high_scores[0]
    exp_lambda = 1 / np.mean(high_scores - high_score)

    hits = []
    # run pwm against the database
    for db_model in db_models_g:
        sim, (start1, end1), (start2, end2), (bg_score, cross_score, pad_score) = \
            model_sim(model, db_model, min_overlap_g)
        if sim < high_score:
            # the score is not in the top scores of the background model
            # this is surely not a significant hit
            continue

        pvalue = highscore_fraction_g * np.exp(- exp_lambda * (sim - high_score))
        evalue = db_size_g * pvalue
        if evalue < evalue_thresh_g:
            hits.append((model['model_id'], db_model['model_id'], sim, evalue, start1, end1, start2, end2,
                         max(bg_score, 0), max(cross_score, 0), pad_score))
    return hits


if __name__ == '__main__':
    main()
