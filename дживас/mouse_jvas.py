import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO

import datetime
import gc

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

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)
    #determineColumns(pd_header)

    all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
               '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    #all_chr = ['X', 'Y']
    #all_chr = ['18']

    total_miss = 0
    total_num = 0
    total_chr = 0

    test_amount_seqs = 140

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

            R1, Y1 = purine_utils.calc_RY(sequence, pos, allele_1)
            R2, Y2 = purine_utils.calc_RY(sequence, pos, allele_2)

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