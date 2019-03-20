#!/usr/bin/env python3
'''
This script is used for comparing a query motif to
motifs in a searching database.
Prerequisite: all motifs in database must be in BaMM format.
The query motif can be either a PWM or a BaMM.
Note: currently it only compares motifs by 0th-order.
'''

import logging
import sys
import argparse
import numpy as np
from multiprocessing import Pool


from utils import calculate_H_model_bg, \
    calculate_H_model, model_sim, update_pwms, \
    update_bamms, parse_meme, parse_bamm_db, parse_bamm_file


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="input MEME-format file with multiple PWMs")
    parser.add_argument('db_path', help="specify the path to database that you want to search for")
    parser.add_argument('output_score_file', default=None, help="name the output file with path")

    parser.add_argument('--query_format', default="PWM",
                        help="declare input format: PWM or BaMM. This needs to be consistent with your input model! "
                             "Default: PWM")
    parser.add_argument('--db_format', default="PWM",
                        help="declare DB format: PWM or BaMM. This needs to be consistent with your input model! "
                             "Default: PWM")
    parser.add_argument('--n_neg_perm', type=int, default=10, help="number of negative permutations. Default: 10")
    parser.add_argument('--highscore_fraction', type=float, default=0.1)
    parser.add_argument('--pvalue_threshold', type=float, default=0.01,
                        help="p-value threshold for output models. Default: 0.01")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--min_overlap', type=int, default=2, help="minimal overlaps between PWMs. Default: 2")
    parser.add_argument('--n_processes', type=int, default=4, help="how many cores are used. Default: 4")


    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    query_file = args.input_file
    target_db_path = args.db_path
    output_score_file = args.output_score_file

    # print out logs
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger_fmt = '%(asctime)s [%(levelname)s]  %(message)s'
    formatter = logging.Formatter(logger_fmt)
    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logger.setLevel(logging.INFO)

    # parse input query file
    query_format = args.query_format
    query_model_set = []
    if query_format == "PWM":
        # parse input meme file
        query_model_set = parse_meme(query_file)['models']
        # pre-compute entropy for all query meme models
        query_models = update_pwms(query_model_set)
    elif query_format == "BaMM":
        # parse input bamm file
        query_model_set = parse_bamm_file(query_file)
        logger.info('There are %s query motif(s) parsed.', len(query_model_set))
        # pre-compute entropy for all query meme models
        query_models = update_bamms(query_model_set)

    else:
        logger.info('Input model file is not recognised. ')

    # parse bamms from the target database
    db_format = args.db_format
    logger.info("DB models are in "+db_format+" format.")

    logger.info('Reading in BaMMs from the target database.')
    target_db = parse_bamm_db(target_db_path)

    # pre-compute entropy for all models in target database
    db_models = update_pwms(target_db)

    db_size = len(db_models)

    logger.info("There are "+str(db_size)+" motif models in the motif DB.")

    # initialize task for paralleling jobs
    def init_workers():
        np.random.seed(args.seed)
        global highscore_fraction_g
        highscore_fraction_g = args.highscore_fraction
        global pvalue_thresh_g
        pvalue_thresh_g = args.pvalue_threshold
        global db_models_g
        db_models_g = db_models
        global db_size_g
        db_size_g = db_size
        global n_neg_perm_g
        n_neg_perm_g = args.n_neg_perm
        global min_overlap_g
        min_overlap_g = args.min_overlap

    logger.info('Queuing %s search job(s)', len(query_models))

    with open(output_score_file, 'w') as out:
        print('model_id', 'db_id', 'p-value', 'e-value', 'sim_score',
              'model_width', sep='\t', file=out)
        with Pool(args.n_processes, initializer=init_workers) as pool:
            jobs = []
            for model in query_models:
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
            shuffle_sim, *_ = model_sim(shuffle_model, db_model, min_overlap_g)
            shuffled_dists.append(shuffle_sim)

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
        sim, *_ = model_sim(model, db_model, min_overlap_g)
        if sim < high_score:
            # the score is not in the top scores of the background model
            # this is surely not a significant hit
            continue

        pvalue = highscore_fraction_g * np.exp(- exp_lambda * (sim - high_score))
        evalue = db_size_g * pvalue

        if pvalue < pvalue_thresh_g:
            hits.append((model['model_id'], db_model['model_id'], pvalue, evalue, sim, db_model['motif_length']) )

    return hits


if __name__ == '__main__':
    main()
