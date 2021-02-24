import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO
import matplotlib.pyplot as plt
import seaborn as sns

import datetime
import gc

from argparse import ArgumentParser

class AlgParams:
    reverse_mode = False
    flexible_endings = False

alg_params = AlgParams()

pairs_p = pd.DataFrame( {'x' : pd.Series( [], dtype=int), 'y' : pd.Series( [], dtype=int)})
pairs_h = pd.DataFrame( {'x' : pd.Series( [], dtype=int), 'y' : pd.Series( [], dtype=int)})

class ColumnsSnp:
    id = ''
    chrom = ''
    pos = ''
    other = ''
    effect = ''

columns = ColumnsSnp()

def reset_df_data():
    global pairs_p
    global pairs_h

    pairs_p = pd.DataFrame({'x': pd.Series([], dtype=int), 'y': pd.Series([], dtype=int)})
    pairs_h = pd.DataFrame({'x': pd.Series([], dtype=int), 'y': pd.Series([], dtype=int)})


def determineColumns(df_snp):
    if 'hm_chrom' in df_snp.columns:
        columns.id = 'hm_rsid'
        columns.chrom = 'hm_chrom'
        columns.pos = 'hm_pos'
        columns.other = 'hm_other_allele'
        columns.effect = 'hm_effect_allele'
    elif 'chromosome' in df_snp.columns:
        columns.id = 'variant_id'
        columns.chrom = 'chromosome'
        columns.pos = 'base_pair_location'
        columns.other = 'other_allele'
        columns.effect = 'effect_allele'
    elif 'CHR' in df_snp.columns:
        columns.id = 'MarkerName'
        columns.chrom = 'CHR'
        columns.pos = 'POS'
        columns.other = 'A2'
        columns.effect = 'A1'

        if ('REF' in df_snp.columns) and ('ALT' in df_snp.columns):
            columns.other = 'REF'
            columns.effect = 'ALT'

    else:
        pass

def RY(nucl):
    if (nucl == 'A') or (nucl == 'G'):
        return 'R'
    else:
        return 'Y'

def calc_hydro_letter(two_letters):
    assert(len(two_letters) == 2)

    if two_letters[0] == 'N' or two_letters[1] == 'N':
        left_chain_start_code = 'N'
        return left_chain_start_code

    if two_letters == 'AA':
        left_chain_start_code = 'k'
    elif two_letters == 'AC':
        left_chain_start_code = 'm'
    elif two_letters == 'AG':
        left_chain_start_code = 'm'
    elif two_letters == 'AT':
        left_chain_start_code = 'f'
    elif two_letters == 'CA':
        left_chain_start_code = 'k'
    elif two_letters == 'CC':
        left_chain_start_code = 'm'
    elif two_letters == 'CG':
        left_chain_start_code = 'p'
    elif two_letters == 'CT':
        left_chain_start_code = 'm'
    elif two_letters == 'GA':
        left_chain_start_code = 'k'
    elif two_letters == 'GC':
        left_chain_start_code = 'm'
    elif two_letters == 'GG':
        left_chain_start_code = 'm'
    elif two_letters == 'GT':
        left_chain_start_code = 'm'
    elif two_letters == 'TA':
        left_chain_start_code = 'p'
    elif two_letters == 'TC':
        left_chain_start_code = 'k'
    elif two_letters == 'TG':
        left_chain_start_code = 'k'
    elif two_letters == 'TT':
        left_chain_start_code = 'k'
    else:
        print(two_letters)
        assert(False)

    return left_chain_start_code

def to_purin(seq):
    purine_seq = []
    purine_seq[:] = seq
    for i in range(0, len(seq) ):
        nucl = seq[i]
        if nucl in ['A', 'G']:
            purine_seq[i] = 'R'
        elif nucl in ['C', 'T']:
            purine_seq[i] = 'Y'
        elif nucl == 'N':
            purine_seq[i] = 'N'
        else:
            assert (False)

    return purine_seq

class RangNucl:
    start = 0
    end = 0
    len = 0

    def __init__(self, start, end):
        self.start = start
        self.end = end
        self.len = self.len()

    def len(self):
        return self.end - self.start

    def __eq__(self, other):
        return self.len == other.len

    def __lt__(self, other):
        return self.len < other.len

def find_purine_chains(purin_seq):
    i = 0
    start = -1

    res = []

    while (i < len(purin_seq)):
        if purin_seq[i] != 'R':
            if start >= 0:
                if (i - start >= 3): # min length
                    res.append( RangNucl(start, i) )
                start = -1
        else:
            if start < 0:
                start = i

        i = i + 1

    return res


def find_hydro_chains(hydro_seq):
    i = 0
    start = -1

    res = []

    while (hydro_seq[i] in ['K', 'M']) and (i < len(hydro_seq) ):
        i = i + 1

    while (i < len(hydro_seq) - 1 ):
        if hydro_seq[i] not in ['K', 'M']:
            if start >= 0:
                if (i - start >= 3): # min length
                    if alg_params.flexible_endings:
                        res.append( RangNucl(start , i ) )
                    else:
                        res.append( RangNucl(start - 1, i + 1) )
                start = -1
        else:
            if start < 0:
                start = i

        i = i + 1

    return res





def to_hydro(seq):
    hydro_seq = []
    hydro_seq[:] = seq[:-1]
    for i in range(0, len(seq) - 1 ):
        hydro_seq[i] = calc_hydro_letter(seq[i: (i + 2) ]).upper()

    return hydro_seq

def set_matrix_value_reverse(a_range, b_range, pairs):
    for i in range(0, a_range.len):
        df_new_val = pd.DataFrame({'x': pd.Series([i + a_range.start, b_range.end - i - 1], dtype=int), 'y': pd.Series([ b_range.end - i - 1, i + a_range.start], dtype=int)})

        pairs = pairs.append( df_new_val)

    return pairs


def set_matrix_value(a_range, b_range, pairs):
    for i in range(0, a_range.len):
        df_new_val = pd.DataFrame({'x': pd.Series([i + a_range.start, i + b_range.start], dtype=int), 'y': pd.Series([i + b_range.start, i + a_range.start], dtype=int)})

        pairs = pairs.append( df_new_val)

    return pairs


def check_seq(seq, a_range, b_range, pairs):
    seq_a = seq[a_range.start:a_range.end]
    seq_b = seq[b_range.start:b_range.end]

    if (seq_a == seq_b):
        pairs = set_matrix_value(a_range, b_range, pairs)
    elif( alg_params.reverse_mode and seq_a == seq_b[::-1] ):
        pairs = set_matrix_value_reverse(a_range, b_range, pairs)




    return pairs





def test():
    data_url = 'http://bit.ly/2cLzoxH'
    gapminder = pd.read_csv(data_url)
    print( np.unique( gapminder['continent'] ))

if __name__ == '__main__':
    #test()
    #exit(0)

    print(datetime.datetime.now())

    config = configparser.ConfigParser()

    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file",  required=False)

    args = parser.parse_args()

    if args.config is not None:
        config.read(args.config)
        print(args.config)
    else:
        config.read('config_t32.ini')

    chrom_dir = config['DEFAULT']['chrom_dir']
    output = config['DEFAULT']['output_file']
    alg_params.reverse_mode = config['DEFAULT'].getboolean('reverse_mode')
    alg_params.flexible_endings = config['DEFAULT'].getboolean('flexible_endings')


    all_chr = ['1']

    for cur_chr in all_chr:

        print(cur_chr)

        chr_filename = 'chr' + cur_chr + '.fa'

        full_filename = os.path.join(chrom_dir, chr_filename)

        if not os.path.exists(full_filename):
            continue

        fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')

        for fasta in fasta_seqs:
            name, sequence = fasta.id, str(fasta.seq)
            break

        num_pages = 50

        common_start = 30000000
        for k in range(0, num_pages ):
            reset_df_data()

            sequence = sequence.upper()
            start_cur = common_start + k * 2000
            end_cur =   common_start + k * 2000 + 2000
            chunk_sequence = sequence[start_cur:end_cur]

            print(start_cur, end_cur)
            #chunk_sequence = 'NnnnnnnnCaaaaaatCgatatatatCaaaaaatGgtgnnnnGtctctcaNnnnggggGtctctcg'
            #chunk_sequence = chunk_sequence.upper()



            purin_seq = to_purin(chunk_sequence)
            hydro_seq = to_hydro(chunk_sequence)

            purine_list = find_purine_chains(purin_seq)
            purine_list = sorted(purine_list)

            hydro_list = find_hydro_chains(hydro_seq)
            hydro_list = sorted(hydro_list)


            for cur_group in range (0, len(purine_list) - 1):
                for i in range(cur_group + 1, len(purine_list) ):
                    if purine_list[cur_group].len == purine_list[i].len:
                        pairs_p = check_seq(chunk_sequence, purine_list[cur_group], purine_list[i], pairs_p)

            for cur_group in range (0, len(hydro_list) - 1):
                for i in range(cur_group + 1, len(hydro_list) ):
                    if hydro_list[cur_group].len == hydro_list[i].len:
                        pairs_h = check_seq(chunk_sequence, hydro_list[cur_group], hydro_list[i], pairs_h)


            pairs_common = pd.merge(pairs_p, pairs_h, how="inner", on=["x", "y"])

            pairs_p['type'] = 'Purine'
            pairs_h['type'] = 'Hydrocode'
            pairs_common['type'] = 'Common'

            pairs_all = pd.concat([pairs_p, pairs_h, pairs_common])

            palette = []

            if pairs_p.shape[0] > 0:
                palette.append('blue')

            if pairs_h.shape[0] > 0:
                palette.append('red')

            if pairs_common.shape[0] > 0:
                palette.append('black')

            fig = plt.figure()
            ax = fig.add_subplot(111)  # The big subplot
            ax.set_aspect('equal', adjustable='box')

            sns.scatterplot(data=pairs_all, x="x", y="y", hue="type",  ax=ax, s=0.15,  palette=palette, legend='full')


            ax.set_xlim(0, len(chunk_sequence))
            ax.set_ylim(0, len(chunk_sequence))

            ax.set_xlabel("x", fontweight='bold')
            ax.set_ylabel("y", fontweight='bold')

            ax.legend(loc='upper left', bbox_to_anchor=(1.04, 1), fontsize = 5)

            description = 'Input File: ' + full_filename + '\n'
            description += 'Start MB: {0} '.format(start_cur / 1000000) + '\n'
            description += 'Chunk : {0}'.format(k) + '\n'

            if alg_params.reverse_mode:
                description += 'Reverse_mode: True \n'

            if alg_params.flexible_endings:
                description += 'Flexible Endings: True \n'


            ax.set_title(description, fontweight='bold', fontsize=4)


            fig.tight_layout()

            output_png = '{}_{}.png'.format(output, k)
            plt.savefig(output_png, dpi=2400, bbox_inches='tight')
            plt.clf()
            plt.close()

    print(datetime.datetime.now())


