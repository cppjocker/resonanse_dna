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

from argparse import ArgumentParser

shoulder_codes = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'GA', 'GC', 'TA']


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

    write_shoulders = True

    test_amount_seqs = 140

    shoulder_codes = pd.read_csv('main3.csv',  delimiter=',')

    for cur_chr in all_chr:

        print(cur_chr)

        df_iter = pd.read_csv(snp_table_file, sep='\t', iterator=True, chunksize=4000, header=0,
                              dtype={columns.chrom: str})

        df_part = pd.DataFrame(columns=pd_header.columns)

        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int, columns.freq : float})

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

        # df_part['arm_AT_l'] = -1
        # df_part['arm_AT_r'] = -1
        # df_part['arm_AT'] = -1

        if write_shoulders:
            df_part['shoulder_l'] = ''
            df_part['shoulder_r'] = ''
            df_part['code_id'] = ''



        if test_only:
            df_part['seq1'] = ''
            df_part['seq2'] = ''

            df_part['seq1_h'] = ''
            df_part['seq2_h'] = ''


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

            km1, seq1, seq1_h = hydro_utils.calc_hydro_by_code(sequence, pos - 1, allele_1, code='hd_07')
            km2, seq2, seq2_h = hydro_utils.calc_hydro_by_code(sequence, pos - 1, allele_2, code='hd_07')


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
                        ( ( allele_1 == first_dinucl[0] and allele_2  == second_dinucl[0] ) or \
                          (allele_2 == first_dinucl[0] and allele_1 == second_dinucl[0]) ) and \
                            (sequence[pos] == first_dinucl[1] ) and ((sequence[pos] == second_dinucl[1] )):

                          shouder1_l, shouder1_r, arm_AT1 = seq_utils.calc_shoulder_metric(sequence, pos - 1, allele_1, code =  allele_1 + sequence[pos], break_code  = break_dinucl, case=1)

                          total1 += 1
                    elif (snp_pos == 2) and \
                        ( ( allele_1 == first_dinucl[1] and allele_2  == second_dinucl[1] ) or \
                          (allele_2 == first_dinucl[1] and allele_1 == second_dinucl[1]) ) and \
                            (sequence[pos - 2] == first_dinucl[0] ) and ((sequence[pos - 2] == second_dinucl[0] )):

                          shouder1_l, shouder1_r, arm_AT1 = seq_utils.calc_shoulder_metric(sequence, pos - 1, allele_1, code = sequence[pos - 2]  + allele_1,  break_code = break_dinucl, case=2)

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




            if freq_eff < 0.5:
                df_part.at[index, 'minA_freq'] = freq_eff
                df_part.at[index, 'maxA_freq'] = 1 - freq_eff

                df_part.at[index, 'maxA_r_len'] = R1
                df_part.at[index, 'maxA_y_len'] = Y1
                df_part.at[index, 'maxA_m'] = max(R1, Y1)

                df_part.at[index, 'minA_r_len'] = R2
                df_part.at[index, 'minA_y_len'] = Y2
                df_part.at[index, 'minA_m'] = max(R2, Y2)

                df_part.at[index, 'maxA_p'] = km1
                df_part.at[index, 'minA_p'] = km2
            else:
                df_part.at[index, 'minA_freq'] = 1 - freq_eff
                df_part.at[index, 'maxA_freq'] = freq_eff

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