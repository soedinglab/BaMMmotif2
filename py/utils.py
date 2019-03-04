#!/usr/bin/env python3

import re
import os
import glob
import operator
import numpy as np
from scipy.special import xlogy
from sklearn.cluster import AffinityPropagation

loge2 = np.log(2)
alpha = 0.8

def calculate_H_model_bg(model, bg):
    H = xlogy(model, model).sum(axis=1) / loge2
    H += xlogy(bg, bg).sum() / loge2
    H *= 0.5
    p_bar = 0.5 * (model + bg)
    H -= xlogy(p_bar, p_bar).sum(axis=1) / loge2
    return H


def calculate_H_model(model):
    H = 0.5 * xlogy(model, model).sum(axis=1) / loge2
    return H


def create_slices(m, n, min_overlap=2):
    # m, n are the lengths of the patterns
    # we demand that n is not longer than m
    assert m >= n

    # obviously it's not possible to overlap in this case
    if n < min_overlap:
        return

    # the shorter pattern can be shifted m - n + 1 inside the longer pattern
    # pad the underhangs in the meanwhile
    # each slice (from left to right) stands for:
    # 1. slice in the overlapped region of m,
    # 2. slice in the overlapped region of n,
    # 3. slice in the underhang region from m on the left,
    # 4. slice in the underhang region from m on the right,
    # 5. slice in the underhang region from n
    for i in range(m - n + 1):
        yield slice(i, i + n), slice(0, n), slice(0, i), slice(i + n + 1, m), slice(0, 0)

    # these are the patterns overlapping the edges with at least min_overlap
    # nucleotides
    for ov in range(min_overlap, n):
        yield slice(0, ov), slice(-ov, None), slice(0, 0), slice(m - ov, m), slice(0, n-ov)
        yield slice(-ov, None), slice(0, ov), slice(0, m - ov), slice(0, 0), slice(ov, n)


def model_sim(model1, model2, min_overlap=2):
    # initialization
    pwm1 = model1['pwm']
    pwm2 = model2['pwm']
    H_model1_bg = model1['H_model_bg']
    H_model2_bg = model2['H_model_bg']
    H_model1 = model1['H_model']
    H_model2 = model2['H_model']

    # my design model2 cannot be longer than model1
    models_switched = False
    if len(pwm1) < len(pwm2):
        pwm1, pwm2 = pwm2, pwm1
        H_model1_bg, H_model2_bg = H_model2_bg, H_model1_bg
        H_model1, H_model2 = H_model2, H_model1
        models_switched = True

    # build model for the reverse complement of the query motif(model1)
    pwm2_rev = pwm2[::-1, ::-1]
    H_model2_rev_bg = H_model2_bg[::-1]
    H_model2_rev = H_model2[::-1]

    scores = []
    contributions = []
    slices = []

    # calculate similarity score between model1 and model2
    for sl1, sl2, uh1l, uh1r, uh2 in create_slices(len(pwm1), len(pwm2), min_overlap):

        # so we want the contributions of the background
        background_score = 0
        background_score += H_model1_bg[sl1].sum()
        background_score += H_model2_bg[sl2].sum()

        # plus the contributions from the padded underhangs
        padding_score = 0
        padding_score += (alpha - 1) * H_model1_bg[uh1l].sum()
        padding_score += (alpha - 1) * H_model1_bg[uh1r].sum()
        padding_score += (alpha - 1) * H_model2_bg[uh2].sum()

        # and the contributions of model1 vs. model2
        cross_score = 0
        cross_score += H_model1[sl1].sum()  # entropy of model1
        cross_score += H_model2[sl2].sum()  # entropy of model2
        p_bar = 0.5 * (pwm1[sl1, :] + pwm2[sl2, :])
        p_bar_entropy = xlogy(p_bar, p_bar) / loge2
        cross_score -= p_bar_entropy.sum()

        # calculate similarity score
        scores.append(alpha * background_score - cross_score + padding_score)

        contributions.append((alpha * background_score, cross_score, padding_score))
        slices.append((sl1, sl2))

    # calculate similarity score between model1_rev and model2
    for sl1, sl2, uh1l, uh1r, uh2 in create_slices(len(pwm1), len(pwm2_rev), min_overlap):
        background_score = 0
        # so we want the contributions of the background
        background_score += H_model1_bg[sl1].sum()
        background_score += H_model2_rev_bg[sl2].sum()

        # plus the contributions from the padded underhangs
        padding_score = 0
        padding_score += (1 - alpha) * H_model1_bg[uh1l].sum()
        padding_score += (1 - alpha) * H_model1_bg[uh1r].sum()
        padding_score += (1 - alpha) * H_model2_rev_bg[uh2].sum()

        cross_score = 0
        # and the contributions of model1 vs. model2
        cross_score += H_model1[sl1].sum()  # entropy of model1
        cross_score += H_model2_rev[sl2].sum()  # entropy of model2

        # cross entropy part
        p_bar = 0.5 * (pwm1[sl1, :] + pwm2_rev[sl2, :])
        p_bar_entropy = xlogy(p_bar, p_bar) / loge2
        cross_score -= p_bar_entropy.sum()

        # calculate similarity score
        scores.append(alpha * background_score - cross_score)

        contributions.append((alpha * background_score, cross_score, padding_score))
        slices.append((sl1, sl2))


    # find the maximum similarity score
    # very neat: https://stackoverflow.com/a/6193521/2272172
    max_index, max_score = max(enumerate(scores), key=operator.itemgetter(1))
    max_slice1, max_slice2 = slices[max_index]

    # gotta love python for that: https://stackoverflow.com/a/13335254/2272172
    start1, end1, _ = max_slice1.indices(len(pwm1))
    start2, end2, _ = max_slice2.indices(len(pwm2))

    contrib = contributions[max_index]

    if models_switched:
        return max_score, (start2 + 1, end2), (start1 + 1, end1), contrib
    else:
        return max_score, (start1 + 1, end1), (start2 + 1, end2), contrib


def update_pwms(models):
    upd_models = []
    for model in models:
        model['pwm'] = np.array(model['pwm'], dtype=float)
        model_length, _ = model['pwm'].shape
        model['bg_freq'] = np.array(model['bg_freq'], dtype=float)
        if 'H_model_bg' not in model or 'H_model' not in model:
            model['H_model_bg'] = calculate_H_model_bg(model['pwm'], model['bg_freq'])
            model['H_model'] = calculate_H_model(model['pwm'])
        else:
            model['H_model_bg'] = np.array(model['H_model_bg'], dtype=float)
            model['H_model'] = np.array(model['H_model'], dtype=float)
        upd_models.append(model)
    return upd_models


def update_bamms(models):
    upd_models = []
    for model in models:
        model['pwm'] = np.array(model['pwm'], dtype=float)
        model_length, _ = model['pwm'].shape
        model['bg_freq'] = np.array(model['bg_freq'], dtype=float)
        if 'H_model_bg' not in model or 'H_model' not in model:
            model['H_model_bg'] = calculate_H_model_bg(model['pwm'], model['bg_freq'])
            model['H_model'] = calculate_H_model(model['pwm'])
        else:
            model['H_model_bg'] = np.array(model['H_model_bg'], dtype=float)
            model['H_model'] = np.array(model['H_model'], dtype=float)
        upd_models.append(model)
    return upd_models


def reduce_pwms(models, min_overlap=2):
    total_models = len(models)
    # calculate similarity score for each model with others
    # build a numpy array of similarity scores
    matrix_sim = np.zeros((total_models, total_models), dtype=float)
    for x_idx, model_1 in enumerate(models, start=0):
        for y_idx, model_2 in enumerate(models, start=0):
            sim_score = model_sim(model_1, model_2, min_overlap)[0]
            matrix_sim[x_idx][y_idx] = sim_score

    # apply affinity propagation to cluster models
    af = AffinityPropagation(affinity='precomputed').fit(matrix_sim)
    af_labels = af.labels_
    #print(af_labels)

    new_models = []

    for rank in range(total_models):
        for idx, model in enumerate(models, start=0):
            if af_labels[idx] == rank:
                new_models.append(models[idx])
                break

    return new_models


def parse_meme(meme_input_file):

    dataset = {}

    with open(meme_input_file) as handle:

        line = handle.readline()

        # check the MEME format version
        if 'MEME version 4' not in line:
            raise ValueError('requires MEME minimal file format version 4')
        else:
            dataset['version'] = line.strip()

        # skip the blank line
        line = handle.readline()

        # read in the ALPHABET info
        dataset['alphabet'] = handle.readline().split()[1]

        # skip the blank line
        line = handle.readline()

        # read in the background letter frequencies
        line = handle.readline()

        # if not given, assign 0.25 to each letter
        bg_freqs = [0.25,0.25,0.25,0.25]
        if line.startswith('Background letter frequencies'):
            bg_toks = handle.readline().split()[1::2]
            bg_freqs = [float(f) for f in bg_toks]
            dataset['bg_freq'] = bg_freqs

        models = []

        for line in handle:
            if line.startswith('MOTIF'):

                model = {}
                model['model_id'] = line.split()[1]
                model['bg_freq'] = bg_freqs

                # read in the information line
                readline = True
                while readline:
                    line = handle.readline()
                    if line.startswith('letter-probability matrix'):
                        readline = False
                info_line = line.rstrip('\n')
                model['info'] = info_line
                width_hit = re.compile('w= (\d+)').search(info_line)
                if not width_hit:
                    raise MalformattedMemeError('could not read motif width')

                # read in the PWM
                pwm_length = int(width_hit.group(1))
                pwm = []

                for i in range(pwm_length):
                    pwm.append([float(p) for p in handle.readline().split()])

                pwm_arr = np.array(pwm, dtype=float)
                bg_arr = np.array(bg_freqs, dtype=float)

                model['pwm'] = pwm
                model['motif_length'] = pwm_length
                model['H_model_bg'] = calculate_H_model_bg(pwm_arr, bg_arr).tolist()
                model['H_model'] = calculate_H_model(pwm_arr).tolist()
                models.append(model)
        dataset['models'] = models
        dataset['bg_freq'] = bg_freqs

    return dataset


def parse_bamm_file(bamm_ifile):

    bamms = []
    bamms.append(parse_single_bamm_from_file(bamm_ifile))
    return bamms


def parse_bamm_db(db_path):

    bamms = []

    bamm_suffix = '.ihbcp'

    for eachdir in os.listdir(db_path):
        if eachdir.endswith(bamm_suffix):
            ifile = os.path.join(db_path, eachdir)
            bamm = parse_single_bamm_from_file(ifile)
            bamms.append(bamm)
        else:
            for ifile in glob.glob(os.path.join(db_path, eachdir, '*', bamm_suffix)):
                bamm = parse_single_bamm_from_file(ifile)
                bamms.append(bamm)

    return bamms


def parse_single_bamm_from_file(bamm_ifile):

    # get the basename of the file
    bn = os.path.splitext(os.path.basename(bamm_ifile))[0]
    dir = os.path.splitext(os.path.dirname(bamm_ifile))[0]
    bg_file = os.path.join(dir, bn.split('_motif_')[0].split('_init')[0] +'.hbcp')

    bamm = {}
    if os.path.isfile(bamm_ifile):
        # read in the model order
        motif_order = 0
        bamm['model_id'] = bn

        with open(bamm_ifile) as bamm_file:
            for line in bamm_file:
                if line[0] != '\n':
                    motif_order = motif_order + 1
                else:
                    break

        # count the motif length
        motif_length = int( sum(1 for line in open(bamm_ifile)) / (motif_order + 1) )

        # read in bamm model
        model = {}
        for k in range(motif_order):
            model[k] = []

        with open(bamm_ifile) as bamm_file:
            for j in range(motif_length):
                for k in range(motif_order):
                    model[k].append( [float(p) for p in bamm_file.readline().split()] )
                # skip the blank line
                bamm_file.readline()

        # convert a bamm array to numpy array
        for k in range(motif_order):
            model[k] = np.array(model[k], dtype=float)

        bamm['model'] = model
        bamm['motif_length'] = motif_length
        bamm['motif_order'] = motif_order-1
        bamm['pwm'] = model[0]

        # get background model frequency
        if os.path.isfile(bg_file):
            with open(bg_file) as bgmodel_file:
                line = bgmodel_file.readline()  # skip the first line for K
                line = bgmodel_file.readline()  # skip the second line for Alpha
                bg_freq = [float(p) for p in bgmodel_file.readline().split()]
            bamm['bg_freq'] = np.array(bg_freq, dtype=float)
        else:
            bamm['bg_freq'] = [0.25,0.25,0.25,0.25]

        # set background model frequency

        bamm['H_model_bg'] = calculate_H_model_bg(bamm['pwm'], bamm['bg_freq']).tolist()
        bamm['H_model'] = calculate_H_model(bamm['pwm']).tolist()
    return bamm


def write_meme(dataset, meme_output_file):
    with open(meme_output_file, "w") as fh:
        print(dataset['version'], file=fh)
        print(file=fh)

        print("ALPHABET= " + dataset['alphabet'], file=fh)
        print(file=fh)

        print("Background letter frequencies", file=fh)

        bg_probs = []
        for idx, nt in enumerate(dataset['alphabet']):
            bg_probs.append(nt)
            bg_probs.append(str(dataset['bg_freq'][idx]))
        print(" ".join(bg_probs), file=fh)
        print(file=fh)

        for model in dataset['models']:
            print("MOTIF {}".format(model['model_id']), file=fh)
            print(model['info'], file=fh)
            pwm = model["pwm"]
            for line in pwm:
                print(" ".join(['{:.4f}'.format(x) for x in line]), file=fh)
            print(file=fh)


def write_bamm(pwm, ofile):
    eps = 1e-16
    with open(ofile, "w") as fh:
        for i in range(len(pwm)):
            print(' '.join(['{:.4e}'.format(x+eps) for x in pwm[i]]) + ' \n', file=fh)


class MalformattedMemeError(ValueError):
    pass
