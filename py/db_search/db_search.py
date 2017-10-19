import argparse
from multiprocessing import Pool
import json
import numpy as np
import logging
import sys
from sklearn.cluster import AffinityPropagation

from utils import calculate_H_model_bg, calculate_H_model, model_sim


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_models')
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

    with open(args.input_models) as model_db:
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
    models_count = len(models)

    def init_workers():
        global db_models_g
        db_models_g = db_models
        global db_size_g
        db_size_g = db_size
        global min_overlap_g
        min_overlap_g = args.min_overlap

    logger.info('Queuing %s search jobs', len(models))

    with open(args.output_file, 'w') as out:
        print('model_id', 'db_id', 'sim_score', 'bg_score', 'cross_score',
              sep='\t', file=out)

        with Pool(args.n_processes, initializer=init_workers) as pool:
            jobs = []
            for model in models:
                job = pool.apply_async(motif_compare, args=(model,))
                jobs.append(job)

            total_jobs = len(jobs)
            motif_motif_matrix = [[]]
            for job_index, job in enumerate(jobs, start=1):
                hits = job.get()
                hits.sort(key=lambda x: x[3])
                motif_motif_score = []
                for hit in hits:
                    print(*hit, sep='\t', file=out)
                    motif_motif_score.append(( hit[2]))
                logger.info('Finished (%s/%s)', job_index, total_jobs)
                motif_motif_matrix.append(motif_motif_score)
                print(job_index)
                print(len(motif_motif_matrix))

            # build a numpy array of similarity scores
            x = np.array(motif_motif_matrix)
            # apply affinity propagation to cluster motifs
            af = AffinityPropagation(affinity='precomputed').fit(x)

            #print(af.labels_)

def motif_compare(model):
    pwm = model['pwm']
    model_id = model['model_id']
    bg_freq = model['bg_freq']
    model_len = len(pwm)
    H_model = model['H_model']
    H_model_bg = model['H_model_bg']
    hits = []
    # run pwm against the database
    for db_model in db_models_g:
        sim, (start1, end1), (start2, end2), (bg_score, cross_score) = model_sim(
            pwm, db_model['pwm'],
            H_model_bg, db_model['H_model_bg'],
            H_model, db_model['H_model'],
            min_overlap=min_overlap_g,
        )

        hits.append((model_id, db_model['model_id'], sim,
                         max(bg_score, 0), max(cross_score, 0)))
    return hits

if __name__ == '__main__':
    main()
