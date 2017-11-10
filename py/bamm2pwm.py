import argparse
from collections import defaultdict
from enum import IntEnum
import numpy as np

IUPAC_ALPHABET_SIZE = 11

#an enum to encode the int representation of the iupac nucleotides
class IUPACNucleotide(IntEnum):
    A = 0
    C = 1
    G = 2
    T = 3
    S = 4
    W = 5
    R = 6
    Y = 7
    M = 8
    K = 9
    N = 10


def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('bamm_file')
    parser.add_argument('meme_file')
    return parser


def main():
    parser = create_parser()
    args = parser.parse_args()

    model = {}
    pwm = []
    motif_length = 0
    line_number = 0
    with open(args.bamm_file) as handle:
        for line in handle:
            line_number += 1
            if len(line.split()) == 0:
                motif_length = motif_length + 1
            if len(line.split()) == 4:
                new_line = np.fromstring(line, dtype=float, sep=' ')
                sum = np.sum(new_line)
                new_line /= sum
                pwm.append(new_line)

    #model_order = int( line_number / motif_length ) - 2
    model['pattern_length'] = motif_length
    model['alphabet'] = "ACGT"
    model['pwm'] = pwm
    model['bg_freqs'] = [0.25,0.25,0.25,0.25]

    # generative a IUPAC string
    representative_iupac_nucleotides = init_representative_map()
    int2char = get_iupac_int2char()
    bg_model = get_bg_model()
    iupac_profiles = init_iupac_profiles(representative_iupac_nucleotides, bg_model)
    IUPAC = get_iupac_string(pwm, iupac_profiles, int2char)
    model['iupac_motif'] = IUPAC

    write_meme(model, args.meme_file)


def write_meme(model, ofile):

    with open(ofile, "w") as fh:
        print("MEME version 4", file=fh)
        print(file=fh)

        print("ALPHABET= " + model['alphabet'], file=fh)
        print(file=fh)

        print("Background letter frequencies", file=fh)
        bg_probs = []
        for idx, nt in enumerate(model["alphabet"]):
            bg_probs.append(nt)
            bg_probs.append(str(model["bg_freqs"][idx]))
        print(" ".join(bg_probs), file=fh)
        print(file=fh)

        print("MOTIF {}".format(model["iupac_motif"]), file=fh)
        print(("letter-probability matrix: alength= {} w= {}")
               .format(len(model['alphabet']), model["pattern_length"]), file=fh)

        for line in model["pwm"]:
            print(" ".join(['{:.4f}'.format(x) for x in line]), file=fh)

        print(file=fh)


#generates the map for the amiguous iupac nucleotides e.g.: N -> A, C, G, T
def init_representative_map():
    representative_iupac_nucleotides = defaultdict(list)

    representative_iupac_nucleotides[IUPACNucleotide.A].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.C].append(IUPACNucleotide.C)
    representative_iupac_nucleotides[IUPACNucleotide.G].append(IUPACNucleotide.G)
    representative_iupac_nucleotides[IUPACNucleotide.T].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.S].append(IUPACNucleotide.C)
    representative_iupac_nucleotides[IUPACNucleotide.S].append(IUPACNucleotide.G)

    representative_iupac_nucleotides[IUPACNucleotide.W].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.W].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.R].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.R].append(IUPACNucleotide.G)

    representative_iupac_nucleotides[IUPACNucleotide.Y].append(IUPACNucleotide.C)
    representative_iupac_nucleotides[IUPACNucleotide.Y].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.M].append(IUPACNucleotide.A)
    representative_iupac_nucleotides[IUPACNucleotide.M].append(IUPACNucleotide.C)

    representative_iupac_nucleotides[IUPACNucleotide.K].append(IUPACNucleotide.G)
    representative_iupac_nucleotides[IUPACNucleotide.K].append(IUPACNucleotide.T)

    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)
    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)
    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)
    representative_iupac_nucleotides[IUPACNucleotide.N].append(IUPACNucleotide.N)

    return representative_iupac_nucleotides


#generates a map to translate the int representation of the iupac nucleotides to chars
def get_iupac_int2char():
    int2char = dict()

    int2char[IUPACNucleotide.A] = 'A'
    int2char[IUPACNucleotide.C] = 'C'
    int2char[IUPACNucleotide.G] = 'G'
    int2char[IUPACNucleotide.T] = 'T'
    int2char[IUPACNucleotide.S] = 'S'
    int2char[IUPACNucleotide.W] = 'W'
    int2char[IUPACNucleotide.R] = 'R'
    int2char[IUPACNucleotide.Y] = 'Y'
    int2char[IUPACNucleotide.M] = 'M'
    int2char[IUPACNucleotide.K] = 'K'
    int2char[IUPACNucleotide.N] = 'N'

    return int2char


# returns a sample bg model; perhaps better to read from an external file?
def get_bg_model():
    bg_model = np.zeros(4)
    bg_model[IUPACNucleotide.A] = 0.2
    bg_model[IUPACNucleotide.C] = 0.3
    bg_model[IUPACNucleotide.G] = 0.3
    bg_model[IUPACNucleotide.T] = 0.2
    return bg_model


# init the profiles for the iupac nucleotides with the given bg_model
def init_iupac_profiles(representative_iupac_nucleotides, bg_model, c=0.2, t=0.7):
    iupac_profiles = np.zeros((IUPAC_ALPHABET_SIZE, 4), np.float)

    for iupac_c in range(IUPAC_ALPHABET_SIZE):
        rep = representative_iupac_nucleotides[iupac_c]
        for a in range(4):
            iupac_profiles[iupac_c][a] += c * bg_model[a]
            for r in rep:
                if a == r:
                    iupac_profiles[iupac_c][a] += t
    return iupac_profiles


#calculates the distance between two profiles; based on the Shannon Entropy?
def calculate_d(profile1, profile2):
    d = 0.0
    for a in range(4):
        d += (profile1[a] - profile2[a]) * (np.log2(profile1[a]) - np.log2(profile2[a]))
    return d


#finds for each profile in the pwm the closest iupac profile
def get_iupac_string(pwm, iupac_profiles, int2char):
    res = []

    pattern_length = len(pwm)
    for i in range(pattern_length):
        min_dist = np.inf
        min_iupac = 0
        for m in range(IUPAC_ALPHABET_SIZE):
            dist = calculate_d(pwm[i], iupac_profiles[m])
            if dist < min_dist:
                min_dist = dist
                min_iupac = m

        res.append(int2char[min_iupac])

    return "".join(res)


if __name__ == '__main__':
    main()
