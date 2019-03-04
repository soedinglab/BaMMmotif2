#!/usr/bin/env python3
'''
This script is used for filtering out motifs from shoot_peng with low occurrences
command for running: python3 cherrypick_peng.py input_file output_file
Prerequisite: input file must be in MEME-format, version 4
(C) Wanwan Ge, 2018-09
'''

import argparse
import shutil
import re


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="input MEME-format file with multiple PWMs")
    parser.add_argument('output_file', help="name the output file with path")
    parser.add_argument('--threshold', type=float, default=0.05,
                        help="Threshold for filtering out low-occurred motifs. Default: 0.05")

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    query_file = args.input_file
    out_file = args.output_file
    threshold = args.threshold

    # parse input query meme file
    model_set = parse_meme(query_file)

    if len(model_set['models']) < 2:
        # no real filtering possible, copy input to output
        shutil.copyfile(query_file, out_file)
        return

    # filter out models with occurrence smaller than certain threshold
    new_models = reduce_pwms(model_set['models'], threshold)

    # update the model set after filtering
    new_model_set = {}
    new_model_set['version']    = model_set['version']
    new_model_set['alphabet']   = model_set['alphabet']
    new_model_set['bg_freq']    = model_set['bg_freq']
    new_model_set['models']     = new_models
    # write out the meme file
    write_meme(new_model_set, out_file)


def parse_meme(meme_input_file):
    dataset = {}

    with open(meme_input_file) as handle:

        line = handle.readline()

        # check the MEME format version
        if 'MEME version 4' not in line:
            raise ValueError('requires MEME minimal file format version 4')
        else:
            dataset['version'] = line.strip()

        models = []

        line = handle.readline()

        if line.startswith('\n'):
            line = handle.readline()

        # read in the ALPHABET info
        if line.startswith('ALPHABET'):
            dataset['alphabet'] = line.split()[1]
            line = handle.readline()
        else:
            dataset['alphabet'] = 'unknown'
            line = handle.readline()

        # skip strands lines
        if line.startswith('strands'):
            line = handle.readline()

        if line.startswith('\n'):
            line = handle.readline()

        if line.startswith('Background letter frequencies'):
            bg_toks = handle.readline().split()[1::2]
            bg_freqs = [float(f) for f in bg_toks]
            dataset['bg_freq'] = bg_freqs
            line = handle.readline()
        else:
            # if not given, assign 0.25 to each letter
            dataset['bg_freq'] = [0.25,0.25,0.25,0.25]
            line = handle.readline()


        if line.startswith('\n'):
            line = handle.readline()

        loop_over = True

        while loop_over:

            if line.startswith('MOTIF'):
                model = {}
                model['model_id'] = line.split()[1]
                #model['bg_freq'] = bg_freqs

                # read in the information line
                readline = True
                while readline:
                    line = handle.readline()
                    if line.startswith('letter-probability matrix'):
                        readline = False
                info_line = line.rstrip('\n')
                model['info'] = info_line
                width_hit = re.compile('w= (\d+)').search(info_line)[1]
                if not width_hit:
                    raise MalformattedMemeError('could not read motif width')

                occur = re.compile('occur= ([-+]?\d*\.\d+|\d+)').search(info_line)[1]

                # read in the PWM
                pwm_length = int(width_hit)
                pwm = []

                for i in range(pwm_length):
                    pwm.append([float(p) for p in handle.readline().split()])
                model['pwm'] = pwm
                model['occur'] = float(occur)
                models.append(model)
                line = handle.readline()
                continue
            elif line.startswith('\n'):
                line = handle.readline()
                continue
            else:
                loop_over = False
        dataset['models'] = models

    return dataset


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


def reduce_pwms(models, threshold = 0.05):
    new_models = []
    for idx, model in enumerate(models, start=0):
        if model['occur'] > threshold:
            new_models.append(model)
    return new_models


if __name__ == '__main__':
    main()