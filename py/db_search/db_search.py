import argparse
from multiprocessing import Pool
import json
import numpy as np
import logging
import sys

from utils import calculate_H_model_bg, calculate_H_model, model_sim


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_models')
    parser.add_argument('model_db')
    parser.add_argument('--n_neg_perm', type=int, default=10)
    parser.add_argument('--highscore_fraction', type=float, default=0.1)
    parser.add_argument('--evalue_threshold', type=float, default=0.1)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--min_overlap', type=int, default=4)
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('output_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger_fmt = '%(asctime)s [%(levelname)s]  %(message)s'
    formatter = logging.Formatter(logger_fmt)
    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    logger.setLevel(logging.INFO)

    def update_models(models):
        upd_models = []
        for model in models:
            model['pwm'] = np.array(model['pwm'], dtype=float)
            model_length, _ = model['pwm'].shape
            if model_length < args.min_overlap:
                logger.warn('model %s with length %s too small for the chosen min_overlap (%s).'
                            ' Please consider lowering the min_overlap threshold.',
                            model['model_id'], model_length, args.min_overlap)
                continue
            model['bg_freq'] = np.array(model['bg_freq'], dtype=float)
            if 'H_model_bg' not in model or 'H_model' not in model:
                model['H_model_bg'] = calculate_H_model_bg(model['pwm'], model['bg_freq'])
                model['H_model'] = calculate_H_model(model['pwm'])
            else:
                model['H_model_bg'] = np.array(model['H_model_bg'], dtype=float)
                model['H_model'] = np.array(model['H_model'], dtype=float)
            upd_models.append(model)
        return upd_models

    with open(args.input_models) as in_models:
        models = update_models(json.load(in_models))

    with open(args.model_db) as model_db:
        db_models = update_models(json.load(model_db))
        db_size = len(db_models)

    rev_models = []
    for model in models:
        rev_model = dict(model)
        rev_model['model_id'] = model['model_id'] + '_rev'
        # reverse complement the pwm
        rev_model['pwm'] = model['pwm'][::-1, ::-1]

        # entropy calculations simply reverse
        rev_model['H_model_bg'] = model['H_model_bg'][::-1]
        rev_model['H_model'] = model['H_model'][::-1]
        rev_models.append(rev_model)

    # intertwine models and rev. complemented models
    models = [model for pair in zip(models, rev_models) for model in pair]

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

    with open(args.output_file, 'w') as out:
        print('model_id', 'db_id', 'simscore', 'e-value',
              'start_query', 'end_query', 'start_hit', 'end_hit', 'bg_score', 'cross_score',
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
    pwm = model['pwm']
    model_id = model['model_id']
    bg_freq = model['bg_freq']
    model_len = len(pwm)
    H_model = model['H_model']
    H_model_bg = model['H_model_bg']

    # step 1: use shuffled pwms to estimate the p-value under the null.
    shuffled_dists = []
    for _ in range(n_neg_perm_g):
        # calculate a locality preserving permutation
        assert model_len > 1
        Z = np.random.normal(0, 1, model_len)
        shuffle_ind = np.argsort(2 * Z - np.arange(model_len))

        shuffle_pwm = pwm[shuffle_ind]
        H_shuffle_bg = calculate_H_model_bg(shuffle_pwm, bg_freq)
        H_shuffle = calculate_H_model(shuffle_pwm)

        for db_model in db_models_g:
            shuf_sim, *_ = model_sim(
                shuffle_pwm, db_model['pwm'],
                H_shuffle_bg, db_model['H_model_bg'],
                H_shuffle, db_model['H_model'],
                min_overlap=min_overlap_g
            )
            shuffled_dists.append(shuf_sim)

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
        sim, (start1, end1), (start2, end2), (bg_score, cross_score) = model_sim(
            pwm, db_model['pwm'],
            H_model_bg, db_model['H_model_bg'],
            H_model, db_model['H_model'],
            min_overlap=min_overlap_g,
        )
        if sim < high_score:
            # the score is not in the top scores of the background model
            # this is surely not a significant hit
            continue

        pvalue = highscore_fraction_g * np.exp(- exp_lambda * (sim - high_score))
        evalue = db_size_g * pvalue
        if evalue < evalue_thresh_g:
            hits.append((model_id, db_model['model_id'], sim, evalue,
                         start1, end1, start2, end2,
                         max(bg_score, 0), max(cross_score, 0)))
    return hits

if __name__ == '__main__':
    main()
