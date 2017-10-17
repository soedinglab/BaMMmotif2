import argparse
import json

import numpy as np
from utils import calculate_H_model, calculate_H_model_bg


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('in_json_db')       # input PWM database
    parser.add_argument('out_json_db')      # output PWM database after pre-computation
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    with open(args.in_json_db) as in_db:
        models = json.load(in_db)

    for model in models:
        if 'H_model_bg' not in model or 'H_model' not in model:
            bg = np.array(model['bg_freq'], dtype=float)
            pwm = np.array(model['pwm'], dtype=float)
            model['H_model_bg'] = calculate_H_model_bg(pwm, bg).tolist()
            model['H_model'] = calculate_H_model(pwm).tolist()

    with open(args.out_json_db, 'w') as out_db:
        json.dump(models, out_db, indent=4, separators=(',', ': '), sort_keys=True)


if __name__ == '__main__':
    main()

