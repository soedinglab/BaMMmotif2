import argparse
import json
import re

import numpy as np
from utils import calculate_H_model, calculate_H_model_bg


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('meme_file')
    parser.add_argument('model_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    with open(args.meme_file) as handle:

        line = handle.readline()
        if line.strip() != 'MEME version 4':
            raise ValueError('requires MEME minimal file format version 4')

        # skip over all optional info
        while line and line != 'Background letter frequencies\n':
            line = handle.readline()

        if line != 'Background letter frequencies\n':
            raise MalformattedMemeError('could not find background frequencies')

        bg_toks = handle.readline().split()[1::2]
        bg_freqs = [float(f) for f in bg_toks]

        # parse pwms
        width_pat = re.compile('w= (\d+)')
        models = []
        for line in handle:
            if line.startswith('MOTIF'):
                model = {}
                model['model_id'] = line.split()[1]
                model['bg_freq'] = bg_freqs

                info_line = handle.readline()
                width_hit = width_pat.search(info_line)
                if not width_hit:
                    raise MalformattedMemeError('could not read motif width')
                pwm_length = int(width_hit.group(1))
                pwm = []
                for i in range(pwm_length):
                    pwm.append([float(p) for p in handle.readline().split()])

                pwm_arr = np.array(pwm, dtype=float)
                bg_arr = np.array(bg_freqs, dtype=float)

                model['pwm'] = pwm
                model['H_model_bg'] = calculate_H_model_bg(pwm_arr, bg_arr).tolist()
                model['H_model'] = calculate_H_model(pwm_arr).tolist()

                models.append(model)

    with open(args.model_file, 'w') as out_db:
        json.dump(models, out_db, indent=4, separators=(',', ': '), sort_keys=True)


class MalformattedMemeError(ValueError):
    pass


if __name__ == '__main__':
    main()
