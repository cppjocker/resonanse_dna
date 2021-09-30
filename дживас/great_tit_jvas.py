import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO

import datetime
import gc

import purine_utils
import hydro_utils
import seq_utils

import re

from argparse import ArgumentParser

class ColumnsSnp:
    id = ''
    chrom = ''
    pos = ''
    other = ''
    effect = ''
    a1_freq = ''
    a2_freq = ''

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
    snp_table_file = 'salmon/salmon_80samples_3-interseq-freq-ancestral.txt'
    output = 'salmon_hydro_hd07_shoulder.tsv'
    chr_filename = 'salmon/GCF_000233375.1_ICSASG_v2_genomic.fna'

    # snp_table_file = 'Great_tit_PRJEB24964/great-tits-vcftools-freq-ancestral-2.txt'
    # output = 'Great_tit_hydro_hd07_shoulders.tsv'
    # chr_filename = 'Great_tit_PRJEB24964/GCF_001522545.3_Parus_major1.1_genomic.fna'


    #snp_table_file = 'wild_rats/Charles-River-freq2.txt'
    # snp_table_file = 'mouse_PMC5020872/Mmc_CAST_freq-ancestral2.txt'
    # output = 'mouse_PMC5020872_hydro.tsv'
    # chr_filename = 'mouse_PMC5020872/genome.fa'

    # snp_table_file = 'wild_rats/Harlan-freq1.txt'
    # output = 'Harlan_hydro.tsv'
    # chr_filename = 'wild_rats/GCA_000001895.4_Rnor_6.0_genomic.fna'


    #snp_table_file = 'bos_taurus/bos_taurus.tsv'
    #output = 'bos_taurus_hydro.tsv'
    #chr_filename = 'bos_taurus/bosTau6.fa'

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)


    total_miss = 0
    total_num = 0
    total_chr = 0

    test_amount_seqs = 140

    write_shoulders = True
    shoulder_codes = pd.read_csv('main3.csv',  delimiter=',')


    chrom_dir = '.'

    full_filename = os.path.join(chrom_dir, chr_filename)

    if not os.path.exists(full_filename):
        assert(False)

    fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')

    #use it to switch between great tit and salmon for example. False = Salmon
    use_chrom_word = False

    for fasta in fasta_seqs:
        name, sequence = fasta.description, str(fasta.seq)
        sequence = sequence.upper()

        print(name)

        if 'random' in name:
            continue

        name = name.replace(',', '')
        words = name.split(' ')

        found = False

        #bos_taurus
        if (len(name) > 3) and (name[0:3] == 'chr'):
            cur_chr = name[3:]
            found = True
        else:
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
        columns.a1_freq = 'A1FREQ'
        columns.a2_freq = 'A2Freq'

        if 'REF_base' in pd_header.columns:
            columns.other = 'REF_base'
        if 'ALT_allele' in pd_header.columns:
            columns.effect = 'ALT_allele'

        if 'ALLELE_FREQ' in pd_header.columns:
            columns.a1_freq = 'ALLELE_FREQ'
            columns.a2_freq = 'ALLELE_FREQ'

        if 'FREQ2' in pd_header.columns:
            columns.a1_freq = 'FREQ1'
            columns.a2_freq = 'FREQ2'


        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int})

            df_part = df_part.append(df_chunk[df_chunk[columns.chrom] == cur_chr])

        print( df_part.shape )

        if( df_part.shape[0] == 0):
            continue

        df_part['freq1'] = df_part[columns.a1_freq]
        df_part['freq2'] = df_part[columns.a2_freq]


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

        if write_shoulders:
            df_part['shoulder_l'] = ''
            df_part['shoulder_r'] = ''
            df_part['code_id'] = ''



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

            km1, seq1, seq1_h = hydro_utils.calc_hydro_by_code(sequence, pos-1, allele_1, 'hd_07')
            km2, seq2, seq2_h = hydro_utils.calc_hydro_by_code(sequence, pos-1, allele_2, 'hd_07')

            if write_shoulders:
                total1 = 0
                total2 = 0

                sloulder_l_all = ''
                sloulder_r_all = ''
                code_i_all = ''

                for code_i in range(shoulder_codes.shape[0]):

                    first_dinucl = shoulder_codes.iloc[code_i, 0]
                    second_dinucl = shoulder_codes.iloc[code_i, 1]
                    break_dinucl = shoulder_codes.iloc[code_i, 2]
                    snp_pos = shoulder_codes.iloc[code_i, 3]

                    if (snp_pos == 1) and \
                            ((allele_1 == first_dinucl[0] and allele_2 == second_dinucl[0]) or \
                             (allele_2 == first_dinucl[0] and allele_1 == second_dinucl[0])) and \
                            (sequence[pos] == first_dinucl[1]) and ((sequence[pos] == second_dinucl[1])):

                        shouder1_l, shouder1_r, arm_AT1 = seq_utils.calc_shoulder_metric(sequence, pos - 1, allele_1,
                                                                                         code=allele_1 + sequence[pos],
                                                                                         break_code=break_dinucl,
                                                                                         case=1)

                        total1 += 1
                    elif (snp_pos == 2) and \
                            ((allele_1 == first_dinucl[1] and allele_2 == second_dinucl[1]) or \
                             (allele_2 == first_dinucl[1] and allele_1 == second_dinucl[1])) and \
                            (sequence[pos - 2] == first_dinucl[0]) and ((sequence[pos - 2] == second_dinucl[0])):

                        shouder1_l, shouder1_r, arm_AT1 = seq_utils.calc_shoulder_metric(sequence, pos - 1, allele_1,
                                                                                         code=sequence[
                                                                                                  pos - 2] + allele_1,
                                                                                         break_code=break_dinucl,
                                                                                         case=2)
                        total2 += 1
                    else:
                        continue

                    sloulder_l_all += '{0},'.format(shouder1_l)
                    sloulder_r_all += '{0},'.format(shouder1_r)
                    code_i_all += '{0},'.format(code_i)

                if total1 + total2 > 0:
                    sloulder_l_all = sloulder_l_all[:-1]
                    sloulder_r_all = sloulder_r_all[:-1]
                    code_i_all = code_i_all[:-1]

                df_part.at[index, 'shoulder_l'] = sloulder_l_all
                df_part.at[index, 'shoulder_r'] = sloulder_r_all
                df_part.at[index, 'code_id'] = code_i_all


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