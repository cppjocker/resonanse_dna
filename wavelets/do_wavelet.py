import pywt
import pandas as pd
import numpy as np

import configparser
from argparse import ArgumentParser
import argparse
import os
import random

import seaborn as sns
import matplotlib.pyplot as plt

import math

from Bio import SeqIO

class AlgParams:
    input_file = ''
    png_file = ''

    do_insert = True
    do_insert_fract = True
    do_randomize = True

    alg_mode = 10
    code_mode = 11 # 0 is 0->1 code, 1 is -1->1 code

    chunk_start = 10
    chunk_end   = 10

class ChromWavelet:
    alg_params = AlgParams()

    def __init__(self, _alg_params):
            self.alg_params = _alg_params

    def get_description(self):
        desc = ''

        params = self.alg_params

        desc = desc + 'Input file path: ' + params.input_file + '. \n'


        if params.do_insert:
            desc += 'Randomized. '

        if params.do_randomize:
            desc += 'Permutted. '

        desc += 'Wavelet spectrogram '

        if params.code_mode == 0:
            desc += 'Purine code [0, 1]. '
        elif params.code_mode == 1:
            desc += 'Purine code [-1, 1]. '

        if params.chunk_start == 0 or params.chunk_end == 0:
            desc += 'Whole range. '
        else:
            desc += '{}-{} Mb range. '.format(int ( params.chunk_start / 1000000 ), int ( params.chunk_end / 1000000) )

        desc = desc + 'Image file: ' + os.path.basename(params.png_file) + '. '


        return desc


    def num_2_nucl(self, num):
        if num == 0:
            return 'A'
        elif num == 1:
            return 'C'
        elif num == 2:
            return 'G'
        elif num == 3:
            return 'T'
        else:
            assert(False)

    def insert_random_next(self, sequence, repeats, start, rand_len):
        num_rand = np.random.binomial(n = 3, p = 0.5, size = rand_len)
        if rand_len == 7:
            num_rand = np.array([1, 0, 0, 0, 1, 0, 0], dtype = np.int32)

        seq_rand = [self.num_2_nucl(x) for x in num_rand]

        sequence = list(sequence)

        seq_rand = np.tile(seq_rand, repeats)
        sequence[start: (start + len(seq_rand) ) ] = seq_rand

        assert(all ( sequence[start: (start + len(seq_rand))] == seq_rand ) )

        sequence = ''.join(sequence)

        return sequence


    def insert_random(self, sequence, repeats):
        insert_place1 =  0
        insert_place2 =  1000000
        insert_place3 = 60000000

        rand_len1 = 190
        rand_len2 = 110
        rand_len3 = 78

        sequence = self.insert_random_next(sequence, repeats, insert_place1, rand_len1)
        sequence = self.insert_random_next(sequence, repeats, insert_place2, rand_len2)
        #sequence = self.insert_random_next(sequence, repeats, insert_place3, rand_len3)

        return sequence


    def make_wavelet(self):
        print(self.alg_params.png_file)

        t = np.arange(1, 1000)
        sig1 = np.sin(2 * math.pi * 1 / 150 * t)
        sig2 = np.sin(2 * math.pi * 1 / 300 * t)

        sig = np.concatenate([sig1, sig2])

        plt.plot(range(0, len(sig)), sig)

        plt.xlabel("Time", fontsize=3)
        plt.ylabel("Amplitude", fontsize=3)

        plt.title("Sin 150 period then 200 period ", fontsize=7)
        plt.savefig("Sin.png", dpi=1200)

        plt.close()

        scales = np.arange(100, 500)

        waveletname = 'morl'
        dt = 1

        [coefficients, frequencies] = pywt.cwt(sig, scales, waveletname, dt)

        idxs_del = np.where(frequencies > 0.5)[0]  # according to Nyquist
        frequencies = np.delete(frequencies, idxs_del)
        coefficients = np.delete(coefficients, idxs_del, 0)

        power = (abs(coefficients)) ** 2
        period = [round(1 / x) for x in frequencies]

        heatmap_pd = pd.DataFrame(
            data=power,  # values
            index=period,
            columns=range(0, len(sig ) ) )

        cs = plt.contourf( range(0, len(sig ) ) , period,  power)
        #sns.set(font_scale=0.5)

        #sns.heatmap(heatmap_pd, cbar_kws={'label': 'Energy'}, cmap='plasma')

        plt.xlabel("Time", fontsize=7)
        plt.ylabel("period (bp)", fontsize=7)

        plt.xticks(fontsize=7)
        plt.yticks(fontsize=7)

        plt.title("Sin example", fontsize=7)

        plt.colorbar()

        plt.savefig(self.alg_params.png_file, dpi=1200)

        plt.close()

        quit()

        fasta_seqs = SeqIO.parse(open(self.alg_params.input_file), 'fasta')

        for fasta in fasta_seqs:
                name, sequence = fasta.id, str(fasta.seq)
                break

        sequence = sequence.upper()

        for i_s in range(0, len(sequence) ):
            if sequence[i_s] != 'N':
                break

        for i_e in range( len(sequence) - 1, 0, -1 ):
            if sequence[i_e] != 'N':
                break

        if self.alg_params.chunk_start != 0 and self.alg_params.chunk_end != 0:
            i_s = max(i_s, self.alg_params.chunk_start)
            i_e = min(i_e, self.alg_params.chunk_end)

        if i_s >= i_e:
            print("Error with chrom range")
            quit(1)

        sequence = sequence[i_s:i_e]

        if self.alg_params.do_insert_fract:
            sequence = self.insert_random2(sequence, 60000)

        if self.alg_params.do_insert:
            sequence = self.insert_random(sequence, 6000)

        if self.alg_params.do_randomize:
            sequence = ''.join(random.sample(sequence, len(sequence)))

        seq = sequence

        seq_np = np.zeros(len(seq), dtype=np.int32)

        if self.alg_params.code_mode == 0:
            for i in range(0, len(seq)):
                if seq[i] in ['C', 'T']:
                    seq_np[i] = 1
                elif seq[i] == 'N':
                    seq_np[i] = -1
                elif seq[i] in ['A', 'G']:
                    seq_np[i] = 0
                else:
                    assert(False)
            idx_gaps = np.where(seq_np == -1)[0]

            #n_gaps = np.count_nonzero(seq_np == -1)
            #n_R = np.count_nonzero(seq_np == 0)
            #n_Y = np.count_nonzero(seq_np == 1)

        elif self.alg_params.code_mode == 1:
            for i in range(0, len(seq)):
                if seq[i] in ['C', 'T']:
                    seq_np[i] = 1
                elif seq[i] in ['A', 'G']:
                    seq_np[i] = -1
                elif seq[i] == 'N':
                    seq_np[i] = 0
                else:
                    assert(False)

            idx_gaps = np.where(seq_np == 0)[0]

        seq_np_cur = seq_np.copy()

        seq_np_cur[idx_gaps] = np.random.binomial(n=1, p=0.5, size=len(idx_gaps))

        if self.alg_params.code_mode == 1:
            seq_np_cur[idx_gaps] = (seq_np_cur[idx_gaps] * 2) - 1

        scales = np.arange(100, 200)
        waveletname = 'morl'
        dt = 1

        desc = self.get_description()

        [coefficients, frequencies] = pywt.cwt(seq_np_cur - np.mean(seq_np_cur), scales, waveletname, dt)

        idxs_del = np.where(frequencies > 0.5)[0] #according to Nyquist
        frequencies = np.delete(frequencies, idxs_del)
        coefficients = np.delete(coefficients, idxs_del, 0)

        power = (abs(coefficients)) ** 2
        period = [ round( 1 / x) for x in frequencies ]
        chrom_poses = [ round(x / 1000000, 1) for x in range(i_s, i_e) ]

        heatmap_pd = pd.DataFrame(
            data=power,  # values
            index=period,
            columns=chrom_poses)

        sns.set(font_scale=0.5)

        sns.heatmap(heatmap_pd, cbar_kws={'label': 'Energy'}, cmap = 'plasma' )

        plt.xlabel("chrom pos (Mb)", fontsize = 7)
        plt.ylabel("period (bp)", fontsize = 7)


        plt.title(desc, fontsize = 7)

        plt.savefig(self.alg_params.png_file, dpi=1200)

        plt.close()




def config2params(config):
    algParams = AlgParams()

    algParams.input_file = os.path.join(config['DEFAULT']['input_dir'],  config['DEFAULT']['input_file'])
    algParams.png_file = os.path.join(config['DEFAULT']['output_dir'],  config['DEFAULT']['png_file'])

    algParams.alg_mode = config['ALG'].getint('alg_mode')
    algParams.code_mode = config['ALG'].getint('code_mode')
    algParams.chunk_start = config['ALG'].getint('chunk_start')
    algParams.chunk_end = config['ALG'].getint('chunk_end')

    algParams.do_insert = config['RANDOM'].getboolean('do_insert')
    algParams.do_insert_fract = config['RANDOM'].getboolean('do_insert_fract')
    algParams.do_randomize = config['RANDOM'].getboolean('do_randomize')

    return algParams

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise NotADirectoryError

if __name__ == '__main__':
    config = configparser.ConfigParser()

    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file", required=False)

    parser.add_argument("-in", "--input",
                        help="chrom Fasta file")

    parser.add_argument("-out", "--output",
                        help="heatmap file", type=argparse.FileType('w'))

    parser.add_argument("-indir", "--input_dir",
                        help="input directory of fasta file", type = dir_path)

    parser.add_argument("-outdir", "--output_dir",
                        help="input directory of fasta file")

    args = parser.parse_args()

    if args.config is not None:
        config.read(args.conf)
    else:
        config.read('config.ini')

    if args.input is not None:
        config['DEFAULT']['input_file'] = args.input

    if args.output is not None:
        config['DEFAULT']['png_file'] = args.out

    if args.input_dir is not None:
        config['DEFAULT']['input_dir'] = args.input_dir

    if args.output_dir is not None:
        config['DEFAULT']['output_dir'] = args.output_dir

    if not os.path.exists(config['DEFAULT']['output_dir'] ) :
        os.mkdir(config['DEFAULT']['output_dir'])

    algParams = config2params(config)
    chromWavelet = ChromWavelet(algParams)
    chromWavelet.make_wavelet()

