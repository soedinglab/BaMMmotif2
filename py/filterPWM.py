'''
This script is for reducing the redundancy of PWMs by filtering similar PWMs out
and only keep one PWM from each cluster.
Prerequisite: input file must be in MEME-format, version 4
'''

import argparse
import logging
import sys
import shutil

from utils import update_models, reduce_pwms, parse_meme, write_meme


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help="input MEME-format file with multiple PWMs")
    parser.add_argument('output_file', help="name the output file with path")
    parser.add_argument('--min_overlap', type=int, default=2, help="minimal overlaps between PWMs. Default: 2")

    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    min_overlap = args.min_overlap
    query_file = args.input_file
    out_file = args.output_file

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

    if len(model_set['models']) < 2:
        # no real filtering possible, copy input to output
        shutil.copyfile(query_file, out_file)
        return

    # pre-compute entropy for all the models
    models = update_models(model_set['models'])

    #filter models using affinity propagation
    new_models = reduce_pwms(models, min_overlap)

    # update the model set after filtering
    new_model_set = {}
    new_model_set['version']    = model_set['version']
    new_model_set['alphabet']   = model_set['alphabet']
    new_model_set['bg_freq']    = model_set['bg_freq']
    new_model_set['models']     = new_models
    # write out the meme file
    write_meme(new_model_set, out_file)


if __name__ == '__main__':
    main()
