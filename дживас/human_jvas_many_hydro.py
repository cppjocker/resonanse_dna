import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO

import datetime
import gc

import purine_utils
import hydro_utils

from argparse import ArgumentParser

code_list = ['hd_07', 'hd_71', 'hd_72', 'hd_73', 'hd_74', 'hd_75', 'hd_76', 'hd_77', 'hd_78', 'hd_79', 'hd_80', 'hd_81', 'hd_82', 'hd_83', 'hd_84', 'hd_85', 'hd_86']

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
        columns.freq = 'hm_effect_allele_frequency'
    elif 'chromosome' in df_snp.columns:
        columns.id = 'variant_id'
        columns.chrom = 'chromosome'
        columns.pos = 'base_pair_location'
        columns.other = 'other_allele'
        columns.effect = 'effect_allele'
        columns.freq = 'effect_allele_frequency'
    elif 'CHR' in df_snp.columns:
        columns.id = 'MarkerName'
        columns.chrom = 'CHR'
        columns.pos = 'POS'
        columns.other = 'A2'
        columns.effect = 'A1'
        columns.freq = 'EAF_A1'

        if ('REF' in df_snp.columns) and ('ALT' in df_snp.columns):
            columns.other = 'REF'
            columns.effect = 'ALT'

    else:
        pass


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

    test_amount_seqs = 140

    for cur_chr in all_chr:

        print(cur_chr)

        df_iter = pd.read_csv(snp_table_file, sep='\t', iterator=True, chunksize=4000, header=0,
                              dtype={columns.chrom: str})

        df_part = pd.DataFrame(columns=pd_header.columns)

        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int})

            df_part = df_part.append(df_chunk[df_chunk[columns.chrom] == cur_chr])

        for code in code_list:
            df_part[code] = ''

        df_part['purine'] = ''
        df_part['remove'] = False

        df_part['seq1'] = ''
        df_part['seq2'] = ''

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

            freq_eff = row[columns.freq]


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

            R1, Y1 = purine_utils.calc_RY(sequence, pos - 1, allele_1)
            R2, Y2 = purine_utils.calc_RY(sequence, pos - 1, allele_2)



            m1 = max(R1, Y1)
            m2 = max(R2, Y2)

            df_part.at[index, 'purine'] = '{0};{1}'.format(m1, m2)

            for code in code_list:
                km1, seq1, seq1_h = hydro_utils.calc_hydro_by_code(sequence, pos - 1, allele_1, code)
                km2, seq2, seq2_h = hydro_utils.calc_hydro_by_code(sequence, pos - 1, allele_2, code)

                df_part.at[index, code] = '{0};{1}'.format(km1, km2)

            seq1 = [ c.lower() for c in sequence[pos-20:pos+20] ]
            seq2 = seq1.copy()

            seq1[19] = allele_1
            seq2[19] = allele_2

            seq1 = "".join(seq1)
            seq2 = "".join(seq2)

            total_num = total_num + 1

            if test_only:
                df_part.at[index, 'seq1'] = seq1
                df_part.at[index, 'seq2'] = seq2

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