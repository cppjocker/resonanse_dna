import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO

import datetime
import gc

from argparse import ArgumentParser


snp_table_file = ''

class ColumnsSnp:
    id = ''
    chrom = ''
    pos = ''
    other = ''
    effect = ''

columns = ColumnsSnp()

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

def calc_hydro(seq, pos, nucl):
    left_start = seq[pos - 1] + nucl
    len_left = 0

    left_chain_start_letter = calc_hydro_letter(left_start)
    next_hydro_letter = left_chain_start_letter

    if next_hydro_letter in {'k', 'm'}:
        len_left = 1

    seq_orig = nucl
    seq_orig = seq[pos - 1].lower() + seq_orig
    seq_returned = next_hydro_letter.upper()

    for i in range(pos - 1, 0, -1):

        next_two_letters = seq[i - 1] + seq[i]
        next_hydro_letter = calc_hydro_letter(next_two_letters)

        seq_returned = next_hydro_letter + seq_returned
        seq_orig = seq[i - 1].lower() + seq_orig

        if next_hydro_letter not in {'k', 'm'}:
            break

        len_left = len_left + 1

    right_start = nucl + seq[pos + 1]
    len_right = 0

    right_chain_start_letter = calc_hydro_letter(right_start)
    next_hydro_letter = right_chain_start_letter

    if next_hydro_letter in {'k', 'm'}:
        len_right = 1


    seq_returned = seq_returned + next_hydro_letter.upper()
    seq_orig =  seq_orig + seq[pos + 1].lower()

    for i in range(pos + 1, len(seq)):
        next_two_letters = seq[i] + seq[i + 1]
        next_hydro_letter = calc_hydro_letter(next_two_letters)

        seq_returned = seq_returned + next_hydro_letter
        seq_orig =  seq_orig + seq[i + 1].lower()

        if next_hydro_letter not in {'k', 'm'}:
            break

        len_right = len_right + 1

    if left_chain_start_letter not in {'k', 'm'} or right_chain_start_letter not in {'k', 'm'}:
        total_len = max(len_left, len_right)
    else:
        total_len = len_left + len_right

    return total_len, seq_orig, seq_returned



def calc_RY(seq, pos, nucl):
    R_chain_len = 0
    Y_chain_len = 0

    ry_nucl = RY(nucl)


    if ry_nucl == 'R':
        R_chain_len = 1
        Y_chain_len = 0
        ry_nucl_compl = 'Y'
    else:
        R_chain_len = 0
        Y_chain_len = 1

        ry_nucl_compl = 'R'

    chain_len = 0

    for i in range(pos - 1, 0, -1):
        if RY(seq[i]) != ry_nucl:
            break

        chain_len+=1

    for i in range(pos + 1, len(seq)):
        if RY(seq[i]) != ry_nucl:
            break

        chain_len+=1

    if ry_nucl == 'R':
        R_chain_len += chain_len
    else:
        Y_chain_len += chain_len

    chain_len_compl_neg = 0
    chain_len_compl_pos = 0

    for i in range(pos - 1, 0, -1):
        if RY(seq[i]) != ry_nucl_compl:
            break

        chain_len_compl_neg+=1

    for i in range(pos + 1, len(seq)):
        if RY(seq[i]) != ry_nucl_compl:
            break

        chain_len_compl_pos+=1

    chain_len_compl = max(chain_len_compl_neg, chain_len_compl_pos)

    if ry_nucl_compl == 'R':
        R_chain_len += chain_len_compl
    else:
        Y_chain_len += chain_len_compl

    return R_chain_len, Y_chain_len

def human_parse():
    print(datetime.datetime.now())

    config = configparser.ConfigParser()

    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file", required=False)

    args = parser.parse_args()

    if args.config is not None:
        config.read(args.config)
        print(args.config)
    else:
        config.read('config.ini')

    chrom_dir = config['DEFAULT']['chrom_dir']
    snp_table_file = config['DEFAULT']['input_file']
    output = config['DEFAULT']['output_file']

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)
    determineColumns(pd_header)

    #    df_snp['seq1'] = ''
    #    df_snp['seq2'] = ''

    #    df_snp['seq1_h'] = ''
    #    df_snp['seq2_h'] = ''

    all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
               '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    #    all_chr = ['18']

    total_miss = 0
    total_num = 0
    total_chr = 0

    for cur_chr in all_chr:

        print(cur_chr)

        df_iter = pd.read_csv(snp_table_file, sep='\t', iterator=True, chunksize=4000, header=0,
                              dtype={columns.chrom: str})

        df_part = pd.DataFrame(columns=pd_header.columns)

        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int})

            df_part = df_part.append(df_chunk[df_chunk[columns.chrom] == cur_chr])

        df_part['r1_len'] = 0
        df_part['r2_len'] = 0

        df_part['y1_len'] = 0
        df_part['y2_len'] = 0

        df_part['m1'] = 0
        df_part['m2'] = 0

        df_part['p1'] = 0
        df_part['p2'] = 0

        df_part['remove'] = False

        chr_filename = 'chr' + cur_chr + '.fa'

        full_filename = os.path.join(chrom_dir, chr_filename)

        if not os.path.exists(full_filename):
            continue

        fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')

        for fasta in fasta_seqs:
            name, sequence = fasta.id, str(fasta.seq)
            break

        sequence = sequence.upper()
        pass

        #       TODO remove this section

        #        text_file = open("seq_test", "w")
        #        text_file.write("%s" % sequence[0:100000])
        #        text_file.close()

        #        quit()

        df_write = []

        for index, row in df_part.iterrows():
            allele_1 = row[columns.other]
            allele_2 = row[columns.effect]

            pos = row[columns.pos]

            nucl = sequence[pos - 1]

            if len(allele_1) > 1 or len(allele_2) > 1:
                df_part.at[index, 'remove'] = True
                continue

            if len(allele_1) > 1:
                print(allele_1)

            if allele_1 != nucl:
                pass

            if not (allele_1 == nucl or allele_2 == nucl):
                print(total_miss)
                total_miss = total_miss + 1
            #                continue

            assert (allele_1 == nucl or allele_2 == nucl)

            R1, Y1 = calc_RY(sequence, pos - 1, allele_1)
            R2, Y2 = calc_RY(sequence, pos - 1, allele_2)

            km1, seq1, seq1_h = calc_hydro(sequence, pos - 1, allele_1)
            km2, seq2, seq2_h = calc_hydro(sequence, pos - 1, allele_2)

            df_part.at[index, 'r1_len'] = R1
            df_part.at[index, 'y1_len'] = Y1

            df_part.at[index, 'm1'] = max(R1, Y1)

            df_part.at[index, 'r2_len'] = R2
            df_part.at[index, 'y2_len'] = Y2

            df_part.at[index, 'm2'] = max(R2, Y2)

            df_part.at[index, 'p1'] = km1
            df_part.at[index, 'p2'] = km2

            total_num = total_num + 1

        #            if total_num >= 200:
        #                break

        #            df_part.at[index, 'seq1'] = seq1
        #            df_part.at[index, 'seq2'] = seq2

        #            df_part.at[index, 'seq1_h'] = seq1_h
        #            df_part.at[index, 'seq2_h'] = seq2_h

        #        if total_num >= 200:
        #            df_part = df_part.iloc[0:200]
        #            break

        print(df_part.shape[0])
        df_part = df_part.drop(df_part[df_part.remove].index)
        print(df_part.shape[0])

        df_part = df_part.drop(columns='remove')

        if total_chr == 0:
            df_part.to_csv(output, sep='\t', index=False)
        else:
            df_part.to_csv(output, sep='\t', index=False, header=False, mode='a')

        total_chr = total_chr + 1

        gc.collect()

    print(datetime.datetime.now())


def split_comma(x):
    #print(x.split(',')[0], x.split(',')[1])
    return x.split(',')[0], x.split(',')[1]


def mouse_parse():
    print(datetime.datetime.now())

    config = configparser.ConfigParser()

    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file", required=False)

    args = parser.parse_args()

    if args.config is not None:
        config.read(args.config)
        print(args.config)
    else:
        config.read('config.ini')

    chrom_dir = config['DEFAULT']['chrom_dir']
    snp_table_file = config['DEFAULT']['input_file']
    output = config['DEFAULT']['output_file']
    test_only = config['DEFAULT'].getboolean('test_only')

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)
    #determineColumns(pd_header)

    #    df_snp['seq1'] = ''
    #    df_snp['seq2'] = ''

    #    df_snp['seq1_h'] = ''
    #    df_snp['seq2_h'] = ''

    all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
               '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    all_chr = ['X', 'Y']
    #all_chr = ['18']

    total_miss = 0
    total_num = 0
    total_chr = 0

    test_amount_seqs = 40

    for cur_chr in all_chr:
        print(cur_chr)

        df_iter = pd.read_csv(snp_table_file, sep='\t', iterator=True, chunksize=4000, header=0,
                              dtype={columns.chrom: str})

        df_part = pd.DataFrame(columns=pd_header.columns)
        columns.pos = 'chromStart'
        columns.chrom = '#chrom'
        columns.other = 'allele1'
        columns.effect = 'allele2'
        columns.ref1 = 'refNCBI'
        columns.ref2 = 'refUCSC'


        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int})

            df_part = df_part.append(df_chunk[df_chunk[columns.chrom] == 'chr' + cur_chr])

        print( df_part.shape )

        if( df_part.shape[0] == 0):
            continue

        df_part['allele1'], df_part['allele2'] = zip( *df_part['alleles'].map(split_comma) )
        df_part['freq1'], df_part['freq2'] = zip( *df_part['alleleFreqs'].map(split_comma) )


        df_part['r1_len'] = 0
        df_part['r2_len'] = 0

        df_part['y1_len'] = 0
        df_part['y2_len'] = 0

        df_part['m1'] = 0
        df_part['m2'] = 0

        df_part['p1'] = 0
        df_part['p2'] = 0

        df_part['remove'] = False

        chr_filename = 'chr' + cur_chr + '.fa'

        full_filename = os.path.join(chrom_dir, chr_filename)

        if not os.path.exists(full_filename):
            continue

        fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')

        for fasta in fasta_seqs:
            name, sequence = fasta.id, str(fasta.seq)
            break

        sequence = sequence.upper()
        pass

        #       TODO remove this section

        #        text_file = open("seq_test", "w")
        #        text_file.write("%s" % sequence[0:100000])
        #        text_file.close()

        #        quit()

        df_write = []

        for index, row in df_part.iterrows():
            allele_1 = row[columns.other]
            allele_2 = row[columns.effect]
            ref1 = row[columns.ref1]
            ref2 = row[columns.ref2]

            pos = row[columns.pos]

            nucl = sequence[pos]

            if (len(allele_1) > 1) or (len(allele_2) > 1) or (ref1 not in {'A', 'C', 'G', 'T'}) or (ref2 not in {'A', 'C', 'G', 'T'} ) \
                    or (allele_1 not in {'A', 'C', 'G', 'T'}) or (allele_2 not in {'A', 'C', 'G', 'T'} ):
                df_part.at[index, 'remove'] = True
                continue

            if len(allele_1) > 1:
                print(allele_1)

            if allele_1 != nucl:
                pass

            #if not (allele_1 == nucl or allele_2 == nucl):
            #    print(total_miss)
            #    total_miss = total_miss + 1
            #                continue

            #print(nucl, allele_1, allele_2, ref1, ref2,  pos)
            assert (allele_1 == nucl or allele_2 == nucl or ref1 == nucl or ref2 == nucl)

            R1, Y1 = calc_RY(sequence, pos, allele_1)
            R2, Y2 = calc_RY(sequence, pos, allele_2)

            km1, seq1, seq1_h = calc_hydro(sequence, pos, allele_1)
            km2, seq2, seq2_h = calc_hydro(sequence, pos, allele_2)

            df_part.at[index, 'r1_len'] = R1
            df_part.at[index, 'y1_len'] = Y1

            df_part.at[index, 'm1'] = max(R1, Y1)

            df_part.at[index, 'r2_len'] = R2
            df_part.at[index, 'y2_len'] = Y2

            df_part.at[index, 'm2'] = max(R2, Y2)

            df_part.at[index, 'p1'] = km1
            df_part.at[index, 'p2'] = km2

            total_num = total_num + 1

            if test_only:
                df_part.at[index, 'seq1'] = seq1
                df_part.at[index, 'seq2'] = seq2

                df_part.at[index, 'seq1_h'] = seq1_h
                df_part.at[index, 'seq2_h'] = seq2_h

                if total_num >= test_amount_seqs:
                    df_part = df_part.iloc[0:test_amount_seqs]
                    break

        print(df_part.shape[0])
        df_part = df_part.drop(df_part[df_part.remove].index)
        print(df_part.shape[0])

        df_part = df_part.drop(columns='remove')

     #   if total_chr == 0:
     #       df_part.to_csv(output, sep='\t', index=False)
     #   else:
        df_part.to_csv(output, sep='\t', index=False, header=False, mode='a')

        total_chr = total_chr + 1

        gc.collect()

        if test_only:
            if total_num >= test_amount_seqs:
                break

    print(datetime.datetime.now())

if __name__ == '__main__':
    #human_parse()
    mouse_parse()



