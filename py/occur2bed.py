#!/usr/bin/env python3
'''
 This script is for converting .occurrence file to .bed file
 The output file contains:
 CHROM START END STRAND
'''


import argparse
import os
import pandas as pd


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('occurence_file')
    parser.add_argument('-o', default=None)
    return parser


def parse_occur(ipath):
    df = pd.read_csv(ipath, sep='\t',header=1, names=['chrom', 'length', 'strand', 'pos', 'pattern', 'p-val', 'e-val'])
    chrom, srange = df['chrom'].str.split(':').str
    _, chrom = chrom.str.split('>').str

    if srange.str.find('-')[1] == -1:
        print("Error: the header line contains no information for generating bed file!")
        exit()
    else:
        sleft, sright = srange.str.split('-').str
        strand = df['strand']
        mstart, _, mend = df['pos'].str.split('.').str
        length = df['length'][1]

        sleft=sleft.astype(int)
        sright=sright.astype(int)
        mstart=mstart.astype(int)
        mend=mend.astype(int)
        start, end = [], []
        for i in range(len(chrom)):
            if strand[i] == '+':
                s_start = sleft[i] + mstart[i] -1
                s_end = sleft[i] + mend[i] - 1
            else:
                s_start = sleft[i] + 2 * length - mend[i] -1
                s_end = sleft[i] + 2 * length - mstart[i] - 1
            start.append(s_start)
            end.append(s_end)

        occurs = pd.DataFrame(columns=['#CHROM', 'START', 'END', 'STRAND', 'Pval'])
        occurs['#CHROM'] = chrom
        occurs['START'] = start
        occurs['END'] = end
        occurs['STRAND'] = strand
        occurs['Pval'] = df['p-val']
        return occurs


def write_bed(occurs, ofile):
    occurs.to_csv(ofile, sep='\t', index=False, float_format='%0.2e')


def main():
    parser = create_parser()
    args = parser.parse_args()

    ipath = args.occurence_file

    if args.o is None:
        dir = os.path.dirname(ipath)
    else:
        dir = args.o

    if not os.path.exists(dir):
        os.makedirs(dir)
    basename = os.path.splitext(os.path.basename(ipath))[0] # note this needs to be changed.
    ofile = os.path.join(dir, basename + ".bed")
    occurs = parse_occur(ipath)

    write_bed(occurs, ofile)


if __name__ == '__main__':
    main()