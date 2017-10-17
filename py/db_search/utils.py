import numpy as np
from scipy.special import xlogy
import operator

loge2 = np.log(2)
alpha = 0.5

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

def create_slices(m, n, min_overlap):

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


def model_sim(model1, model2, H_model1_bg, H_model2_bg, H_model1, H_model2, min_overlap=2):

    models_switched = False

    # my design model2 cannot be longer than model1
    if len(model1) < len(model2):
        model1, model2 = model2, model1
        H_model1_bg, H_model2_bg = H_model2_bg, H_model1_bg
        H_model1, H_model2 = H_model2, H_model1
        models_switched = True


    scores = []
    contributions = []
    slices = []

    for sl1, sl2, uh1l, uh1r, uh2 in create_slices(len(model1), len(model2), min_overlap):
        background_score = 0
        # so we want the contributions of the background
        background_score += H_model1_bg[sl1].sum()
        background_score += H_model2_bg[sl2].sum()

        # plus the contributions from the padded underhangs
        padding_score = 0
        padding_score += (1 - alpha) * H_model1_bg[uh1l].sum()
        padding_score += (1 - alpha) * H_model1_bg[uh1r].sum()
        padding_score += (1 - alpha) * H_model2_bg[uh2].sum()

        cross_score = 0
        # and the contributions of model1 vs. model2
        cross_score += H_model1[sl1].sum()  # entropy of model1
        cross_score += H_model2[sl2].sum()  # entropy of model2

        # cross entropy part
        p_bar = alpha * (model1[sl1, :] + model2[sl2, :])
        p_bar_entropy = xlogy(p_bar, p_bar) / loge2
        cross_score -= p_bar_entropy.sum()

        scores.append(background_score - cross_score)
        contributions.append((background_score, cross_score))
        slices.append((sl1, sl2))

    # very neat: https://stackoverflow.com/a/6193521/2272172
    max_index, max_score = max(enumerate(scores), key=operator.itemgetter(1))
    max_slice1, max_slice2 = slices[max_index]

    # gotta love python for that: https://stackoverflow.com/a/13335254/2272172
    start1, end1, _ = max_slice1.indices(len(model1))
    start2, end2, _ = max_slice2.indices(len(model2))

    contrib = contributions[max_index]

    if models_switched:
        return max_score, (start2 + 1, end2), (start1 + 1, end1), contrib
    else:
        return max_score, (start1 + 1, end1), (start2 + 1, end2), contrib
