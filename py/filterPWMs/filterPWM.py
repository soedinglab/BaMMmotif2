import argparse
import json
import numpy as np
import sys

from utils import parse_meme, write_meme, update_models, filter_pwms

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', help='input meme file with pwms')
    parser.add_argument('--min_overlap', type=int, default=4, help='minimal overlap in nt')
    parser.add_argument('--n_processes', type=int)
    parser.add_argument('output_file', help='output file for storing representative pwms')
    args = parser.parse_args()

    if args.output_file is None:
        print("Warning: you did not define an output file (options -o ). Stopping here.", file = sys.stderr)
        sys.exit(1)

    # parse input meme file
    model_set = parse_meme(args.input_file)

    # pre-compute for all the models
    models = update_models(model_set['models'])

    # filter models using affinity propagation
    new_models = filter_pwms(models)

    # update the model set after filtering
    new_model_set = {}
    new_model_set['version'] = model_set['version']
    new_model_set['alphabet'] = model_set['alphabet']
    new_model_set['bg_freq'] = model_set['bg_freq']
    new_model_set['models'] = new_models

    # write out the meme file
    write_meme(new_model_set, args.output_file)


if __name__ == '__main__':
    main()
