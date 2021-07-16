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

code_list = ['hd_07', 'hd_55', 'hd_56', 'hd_57', 'hd_58', 'hd_59', 'hd_60', 'hd_61', 'hd_62', 'hd_63', 'hd_64', 'hd_65',
             'hd_66', 'hd_67', 'hd_68', 'hd_69', 'hd_70']


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


    test_only = False
    # snp_table_file = 'salmon/salmon_80samples_3-interseq-freq-ancestral.txt'
    # output = 'salmon_hydro_55-70.tsv'
    # chr_filename = 'salmon/GCF_000233375.1_ICSASG_v2_genomic.fna'

    snp_table_file = 'Great_tit_PRJEB24964/great-tits-vcftools-freq-ancestral-2.txt'
    output = 'Great_tit_hydro_55-70.tsv'
    chr_filename = 'Great_tit_PRJEB24964/GCF_001522545.3_Parus_major1.1_genomic.fna'


    #snp_table_file = 'wild_rats/Charles-River-freq2.txt'
    # snp_table_file = 'wild_rats/Harlan-freq1.txt'
    # output = 'Harlan_hydro.tsv'
    # chr_filename = 'wild_rats/GCA_000001895.4_Rnor_6.0_genomic.fna'


    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)


    total_miss = 0
    total_num = 0
    total_chr = 0

    test_amount_seqs = 140

    chrom_dir = '.'

    full_filename = os.path.join(chrom_dir, chr_filename)

    if not os.path.exists(full_filename):
        assert(False)

    fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')

    #use it to switch between great tit and salmon for example
    # True -> GreatTit. False->Salmon
    use_chrom_word = True

    for fasta in fasta_seqs:
        name, sequence = fasta.description, str(fasta.seq)
        sequence = sequence.upper()

        print(name)

        if 'random' in name:
            continue

        name = name.replace(',', '')
        words = name.split(' ')


        found = False
        for i in range(0, len(words)):
            if(words[i] == 'chromosome'):
                if use_chrom_word:
                    cur_chr = words[i + 1]
                else:
                    cur_chr = words[0]
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

        for code in code_list:
            df_part[code] = ''
            if test_only:
                df_part[code + "_seq1"] = ''
                df_part[code + "_seq2"] = ''

        df_part['purine'] = ''
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
                df_part.at[index, 'remove'] = True

                print(total_miss)
                total_miss = total_miss + 1
                continue

            assert (allele_1 == nucl or allele_2 == nucl )

            R1, Y1 = purine_utils.calc_RY(sequence, pos-1, allele_1)
            R2, Y2 = purine_utils.calc_RY(sequence, pos-1, allele_2)

            m1 = max(R1, Y1)
            m2 = max(R2, Y2)

            df_part.at[index, 'purine'] = '{0};{1}'.format(m1, m2)

            for code in code_list:
                km1, seq1, seq1_h = hydro_utils.calc_hydro_by_code(sequence, pos - 1, allele_1, code)
                km2, seq2, seq2_h = hydro_utils.calc_hydro_by_code(sequence, pos - 1, allele_2, code)

                df_part.at[index, code] = '{0};{1}'.format(km1, km2)

                if test_only:
                    df_part.at[index, 'seq1'] = seq1
                    df_part.at[index, 'seq2'] = seq2

                    df_part.at[index, code + "_seq1"] = seq1_h
                    df_part.at[index, code + "_seq2"] = seq2_h


            total_num = total_num + 1

            if test_only:
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