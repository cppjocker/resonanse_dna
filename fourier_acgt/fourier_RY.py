from Bio import SeqIO
import numpy as np
import pandas as pd
import math
from scipy import signal
import seaborn as sns
import matplotlib.pyplot as plt
import random

import ctypes
import glob

import configparser
from argparse import ArgumentParser
import argparse
import os


class AlgParams:
    input_file = ''
    heatmap_file = ''

    do_insert = True
    do_insert_fract = True
    do_randomize = True

    alg_mode = 10
    code_mode = 11 # 0 is 0->1 code, 1 is -1->1 code

    chunk_start = 10
    chunk_end   = 10


class ChromFourier:
    alg_params = AlgParams()

    def __init__(self, _alg_params):
        self.alg_params = _alg_params

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


    def make_periodogramm(self, seq_np, size2, idx_gaps):
        size2_half = int(size2 / 2) + 1

        n_exps = 1

        n_gaps = len(idx_gaps)

        if n_gaps > 0:
            n_exps = 200

        psd = np.zeros(size2_half)

        for i in range(0, n_exps):
            seq_np_cur = np.copy(seq_np)
            seq_np_cur[idx_gaps] = np.random.binomial(n=1, p=0.5, size=n_gaps)

            if self.alg_params.code_mode == 1:
                seq_np_cur[idx_gaps] = (seq_np_cur[idx_gaps] * 2) - 1

            seq_np_cur = seq_np_cur.astype(np.float64)

            f, pxx = signal.periodogram(x=seq_np_cur, nfft=size2)
            psd += pxx

        psd /= n_exps

        return psd

    def make_correlation(self, a, b):
        if (a == -1) or (b == -1):
            return np.nan

        if a == b:
            return 1.0

        return 0.0

    def make_corr_func_fast(self, seq_np, size2):
        len_corr_func = int(len(seq_np) / 2)

        corr_func = np.zeros(len_corr_func)

        for i in range(0, len_corr_func):
            end_pos = len(seq_np) - i

            seq_orig = seq_np[0:end_pos].copy()
            seq_shift = seq_np[i:(end_pos+i)].copy()

            vect = np.vectorize(self.make_correlation)
            res_corr = vect(seq_orig, seq_shift)

            total_sum = np.nansum(res_corr)
            total_N = np.count_nonzero(~np.isnan(res_corr))

            if total_N > 0.000001:
                corr_func[i] = total_sum / total_N

        psd_cmplx = np.fft.fft(corr_func - np.mean(corr_func), size2)
        psd = [np.real(x) for x in psd_cmplx]

        size2_half = int(len(psd) / 2 ) + 1

        psd = psd[0:size2_half]

        return psd


    def make_corr_func(self, seq_np, size2):
        len_corr_func = int(len(seq_np) / 2)

        corr_func = np.zeros(len_corr_func)

        for i in range(0, len_corr_func):
            total_sum = 0
            total_N = 0

            for k in range(0, len_corr_func):
                if k + i >= len(seq_np):
                    break

                if (seq_np[k] == -1) or (seq_np[k + i] == -1):   # -1 is gap indicator
                    continue

                total_N += 1

                if seq_np[k] == seq_np[k + i]:
                    total_sum += 1

            if total_N > 0:
                corr_func[i] = total_sum / total_N

        psd_cmplx = np.fft.fft(corr_func - np.mean(corr_func), size2)
        psd = [np.real(x) for x in psd_cmplx]

        size2_half = int(len(psd) / 2 ) + 1

        psd = psd[0:size2_half]

        return psd

    def make_corr_func_cpp(self, seq_np, size2):
        libfile = glob.glob('./corr_func*.so')[0]
        mylib = ctypes.CDLL(libfile)

        mylib.calc_corr_dna.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_int,
                                np.ctypeslib.ndpointer(dtype=np.float64), ctypes.c_int]

        len_corr_func = int(len(seq_np) / 2)
        corr_func = np.zeros(len_corr_func)

        print(seq_np.dtype)
        print(corr_func.dtype)

        seq_np_float = seq_np.astype(np.float64)

        mylib.calc_corr_dna(seq_np_float, len(seq_np_float), corr_func, len_corr_func)

        psd_cmplx = np.fft.fft(corr_func - np.mean(corr_func), size2)
        psd = [np.real(x) for x in psd_cmplx]

        size2_half = int(len(psd) / 2 ) + 1

        psd = psd[0:size2_half]

        return psd


    def insert_random2(self, sequence, repeats):
        insert_place1 = 38000000
        rand_len1 = 7

        sequence = self.insert_random_next(sequence, repeats, insert_place1, rand_len1)

        return sequence

    def insert_random(self, sequence, repeats):
        insert_place1 =  5000000
        insert_place2 = 25000000
        insert_place3 = 60000000

        rand_len1 = 76
        rand_len2 = 77
        rand_len3 = 78

        sequence = self.insert_random_next(sequence, repeats, insert_place1, rand_len1)
        sequence = self.insert_random_next(sequence, repeats, insert_place2, rand_len2)
        sequence = self.insert_random_next(sequence, repeats, insert_place3, rand_len3)

        return sequence

    #alg_mode:
    #
    #0 - periodogram
    #1 - correlational func
    #
    #

    def make_heatmap(self):

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

        bin = 30000

        bins = math.floor ( len(sequence) / bin)
        size2 = 2 ** math.ceil(math.log(bin) / math.log(2))
        size2_half = int(size2 / 2) + 1

        heatmap_mat = np.zeros( (size2_half, bins ) )

        for k in range(0, bins):

            start_seq = k * 30000
            end_seq = k * 30000 + 30000

            seq = sequence[start_seq:end_seq]


            seq_np = np.zeros(len(seq), dtype = np.int32)

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



            print(len(seq))

            print( (start_seq, end_seq) )

            if self.alg_params.alg_mode == 0:
                psd = self.make_periodogramm(seq_np, size2, idx_gaps)
            elif self.alg_params.alg_mode == 1:
                psd = self.make_corr_func(seq_np, size2)
            elif self.alg_params.alg_mode == 2:
                psd = self.make_corr_func_cpp(seq_np, size2)
                #psd = make_corr_func_fast(seq_np, size2)


            heatmap_mat[:, k] = psd
        #plt.semilogy(f, pxx)
        #plt.xlabel('frequency [Hz]')
        #plt.ylabel('PSD [V**2/Hz]')
        #plt.show()

        #np.savetxt(self.alg_params.chrom_start_file, np.array([i_s]) )

        np.savetxt(self.alg_params.heatmap_file, heatmap_mat, delimiter=";")

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise NotADirectoryError

def config2params(config):
    algParams = AlgParams()

    algParams.input_file = os.path.join(config['DEFAULT']['input_dir'],  config['DEFAULT']['input_file'])
    algParams.heatmap_file = os.path.join(config['DEFAULT']['output_dir'],  config['DEFAULT']['heatmap_file'])

    algParams.alg_mode = config['ALG'].getint('alg_mode')
    algParams.code_mode = config['ALG'].getint('code_mode')
    algParams.chunk_start = config['ALG'].getint('chunk_start')
    algParams.chunk_end = config['ALG'].getint('chunk_end')

    algParams.do_insert = config['RANDOM'].getboolean('do_insert')
    algParams.do_insert_fract = config['RANDOM'].getboolean('do_insert_fract')
    algParams.do_randomize = config['RANDOM'].getboolean('do_randomize')

    return algParams

if __name__ == '__main__':
    config = configparser.ConfigParser()

    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file", type=argparse.FileType('r'), required=False)

    parser.add_argument("-in", "--input",
                        help="chrom Fasta file", type=argparse.FileType('r'))

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
        config['DEFAULT']['heatmap_file'] = args.out

    if args.input_dir is not None:
        config['DEFAULT']['input_dir'] = args.input_dir

    if args.output_dir is not None:
        config['DEFAULT']['output_dir'] = args.output_dir

    if not os.path.exists(config['DEFAULT']['output_dir'] ) :
        os.mkdir(config['DEFAULT']['output_dir'])

    algParams = config2params(config)
    fourierProc = ChromFourier(algParams)
    fourierProc.make_heatmap()

    quit()


    size2 = 2 ** math.ceil(math.log(30000) / math.log(2))
    size2_half = int(size2 / 2) + 1


    #"heatmap.csv"

    heatmap_mat = np.loadtxt(heatmap_file, delimiter=";")
    chr_offset  = np.loadtxt(chrom_start_file).tolist()

    print("Loaded")

    f_ranges = [(0.0, 0.002), (0.002, 0.004), (0.004, 0.007), (0.007, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04), (0.04, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)  ]

    percentile = 83
    n_chr_ranges = 1
    total_chr_chunks = heatmap_mat.shape[1]

    visual_chunk = math.floor( total_chr_chunks / n_chr_ranges )

    for j in range(0, n_chr_ranges + 1):
        chr_start = j * visual_chunk
        chr_end = j * visual_chunk + visual_chunk

        if j == n_chr_ranges:
            chr_end = total_chr_chunks

        if chr_start >= total_chr_chunks:
            break

        for k in range(0, len(f_ranges ) ):
            f_range = f_ranges[k]
            f_start = f_range[0]
            f_end = f_range[1]

            start_q = int(f_start * size2)
            start_q = max(1, start_q)
            end_q = int(f_end * size2)

            heatmap_mat_chunk = heatmap_mat[start_q:end_q, chr_start:chr_end]

            #change matrix by percentile logic

            for i_hmap in range(0, heatmap_mat_chunk.shape[0]):
                perc_xx = np.percentile(heatmap_mat_chunk[i_hmap, :], percentile)
                idxs_more_perc_xx = np.where(heatmap_mat_chunk[i_hmap, :] > perc_xx)[0]
                heatmap_mat_chunk[i_hmap, idxs_more_perc_xx] = perc_xx
                pass

            print(heatmap_mat_chunk.shape)
            print(heatmap_mat_chunk.dtype)

            font_scale = 0.5

#            if(heatmap_mat_chunk.shape[0] < 40):
#                font_scale = 0.2

            round_val = 1


            tmp_indexes = list( range(start_q, end_q) )
            indexes = [round(size2 / x, round_val ) for x in tmp_indexes]

            tmp_indexes_chr = list( range(chr_start, chr_end) )
            indexes_chr = [ round( (30000 * x + chr_offset ) / 1000000, 1) for x in tmp_indexes_chr]

            print(indexes)

            df = pd.DataFrame(data=heatmap_mat_chunk, index=indexes, columns=indexes_chr)


            sns.set(font_scale=font_scale)
            sns.heatmap(df, cbar_kws={'label': 'Amplitude in FFT'}, cmap = 'plasma')

            #sns.color_palette("turbo", as_cmap=True)


            plt.xlabel("chrom pos (Mb)")
            plt.ylabel("period (bp)")

            plt.title("Project - Hotmap1. Program fourier_RY.py version 1.0 (heatmap with 30 kbp bin) ." + chrom_name )

            plt.savefig("heatmap.period_{}-{}.chr_{}.png".format(indexes[-1], indexes[1], j), dpi = 600 )

            plt.close()

    pass