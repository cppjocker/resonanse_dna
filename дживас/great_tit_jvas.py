import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO

import datetime
import gc

import purine_utils
import hydro_utils

import re

from argparse import ArgumentParser

class ColumnsSnp:
    id = ''
    chrom = ''
    pos = ''
    other = ''
    effect = ''

columns = ColumnsSnp()


def split_comma(x):
    #print(x.split(',')[0], x.split(',')[1])
    return x.split(',')[0], x.split(',')[1]



def parse():
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


    test_only = True
    snp_table_file = 'Great_tit_PRJEB24964/great-tits-vcftools-freq-ancestral-2.txt'
    output = 'Great_tit_hydro_test.tsv'

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)


    total_miss = 0
    total_num = 0
    total_chr = 0

    test_amount_seqs = 140

    chrom_dir = '.'
    chr_filename = 'Great_tit_PRJEB24964/GCF_001522545.3_Parus_major1.1_genomic.fna'

    full_filename = os.path.join(chrom_dir, chr_filename)

    if not os.path.exists(full_filename):
        assert(False)

    fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')


    for fasta in fasta_seqs:
        name, sequence = fasta.description, str(fasta.seq)
        sequence = sequence.upper()

        print(name)

        name = name.replace(',', '')
        words = name.split(' ')


        found = False
        for i in range(0, len(words)):
            if(words[i] == 'chromosome'):
                cur_chr = words[i + 1]
                found = True
                break

        if not found:
            print('NOT FOUND!')
            continue

        print(cur_chr)

        df_iter = pd.read_csv(snp_table_file, sep='\t', iterator=True, chunksize=4000, header=0,
                              dtype={columns.chrom: str})

        df_part = pd.DataFrame(columns=pd_header.columns)
        columns.pos = 'POS'
        columns.chrom = 'CHROM'
        columns.other = 'A1'
        columns.effect = 'A2'


        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int})

            df_part = df_part.append(df_chunk[df_chunk[columns.chrom] == cur_chr])

        print( df_part.shape )

        if( df_part.shape[0] == 0):
            continue

        df_part['freq1'] = df_part['A1FREQ']
        df_part['freq2'] = df_part['A2Freq']


        df_part['r1_len'] = 0
        df_part['r2_len'] = 0

        df_part['y1_len'] = 0
        df_part['y2_len'] = 0

        df_part['m1'] = 0
        df_part['m2'] = 0

        df_part['p1'] = 0
        df_part['p2'] = 0

        df_part['minA_freq'] = 0.0
        df_part['maxA_freq'] = 0.0


        df_part['minA_r_len'] = 0
        df_part['maxA_r_len'] = 0

        df_part['minA_y_len'] = 0
        df_part['maxA_y_len'] = 0

        df_part['minA_m'] = 0
        df_part['maxA_m'] = 0

        df_part['minA_p'] = 0
        df_part['maxA_p'] = 0

        df_part['remove'] = False


        for index, row in df_part.iterrows():
            allele_1 = row[columns.other]
            allele_2 = row[columns.effect]

            pos = row[columns.pos]

            nucl = sequence[pos-1]

            if (len(allele_1) > 1) or (len(allele_2) > 1) or (allele_1 not in {'A', 'C', 'G', 'T'}) or (allele_2 not in {'A', 'C', 'G', 'T'} ):
                df_part.at[index, 'remove'] = True
                continue



            if not (allele_1 == nucl or allele_2 == nucl):
                print(allele_1, allele_2, nucl, pos)
            #    print(total_miss)
            #    total_miss = total_miss + 1
            #                continue

            assert (allele_1 == nucl or allele_2 == nucl )

            R1, Y1 = purine_utils.calc_RY(sequence, pos-1, allele_1)
            R2, Y2 = purine_utils.calc_RY(sequence, pos-1, allele_2)

            km1, seq1, seq1_h = hydro_utils.calc_hydro(sequence, pos-1, allele_1)
            km2, seq2, seq2_h = hydro_utils.calc_hydro(sequence, pos-1, allele_2)

            df_part.at[index, 'r1_len'] = R1
            df_part.at[index, 'y1_len'] = Y1

            df_part.at[index, 'm1'] = max(R1, Y1)

            df_part.at[index, 'r2_len'] = R2
            df_part.at[index, 'y2_len'] = Y2

            df_part.at[index, 'm2'] = max(R2, Y2)

            df_part.at[index, 'p1'] = km1
            df_part.at[index, 'p2'] = km2

            if row['freq1'] > row['freq2']:
                df_part.at[index, 'minA_freq'] = row['freq2']
                df_part.at[index, 'maxA_freq'] = row['freq1']

                df_part.at[index, 'maxA_r_len'] = R1
                df_part.at[index, 'maxA_y_len'] = Y1

                df_part.at[index, 'maxA_m'] = max(R1, Y1)

                df_part.at[index, 'minA_r_len'] = R2
                df_part.at[index, 'minA_y_len'] = Y2

                df_part.at[index, 'minA_m'] = max(R2, Y2)

                df_part.at[index, 'maxA_p'] = km1
                df_part.at[index, 'minA_p'] = km2
            else:
                df_part.at[index, 'minA_freq'] = row['freq1']
                df_part.at[index, 'maxA_freq'] = row['freq2']

                df_part.at[index, 'maxA_r_len'] = R2
                df_part.at[index, 'maxA_y_len'] = Y2

                df_part.at[index, 'maxA_m'] = max(R2, Y2)

                df_part.at[index, 'minA_r_len'] = R1
                df_part.at[index, 'minA_y_len'] = Y1

                df_part.at[index, 'minA_m'] = max(R1, Y1)

                df_part.at[index, 'maxA_p'] = km2
                df_part.at[index, 'minA_p'] = km1


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

        if total_chr == 0:
            df_part.to_csv(output, sep='\t', index=False)
        else:
            df_part.to_csv(output, sep='\t', index=False, header=False, mode='a')

        total_chr = total_chr + 1

        gc.collect()

        if test_only:
            if total_num >= test_amount_seqs:
                break


    print(datetime.datetime.now())