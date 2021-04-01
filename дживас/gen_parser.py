from gtfparse import read_gtf
import pandas as pd
import numpy as np
import configparser
import gc

import datetime
import jvas_utils as ju

from argparse import ArgumentParser


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


if __name__ == '__main__':
    print(datetime.datetime.now())


    config = configparser.ConfigParser()


    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file",  required=False)

    args = parser.parse_args()

    if args.config is not None:
        config.read(args.config)
        print(args.config)
    else:
        config.read('config_gen.ini')

    chrom_dir = config['DEFAULT']['chrom_dir']
    wig_dir = config['DEFAULT']['wig_dir']

    snp_table_file = config['DEFAULT']['input_file']
    output = config['DEFAULT']['output_file']
    gtf_file = config['DEFAULT']['gtf_file']

    df_gtf = read_gtf(gtf_file)

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)
    determineColumns(pd_header)


    all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
               '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    #all_chr = ['18']

    total_chr = 0

    for cur_chr in all_chr:
        print(cur_chr)

        gtf_chr = 'chr' + cur_chr

        df_iter = pd.read_csv(snp_table_file, sep = '\t', iterator=True, chunksize=4000, header=0, dtype = {columns.chrom: str} )

        df_part = pd.DataFrame(columns=pd_header.columns)

        for df_chunk in df_iter:
            df_chunk = df_chunk.dropna(subset=[columns.pos])

            df_chunk = df_chunk.astype({columns.chrom: str, columns.pos: int})

            df_part = df_part.append( df_chunk[df_chunk[columns.chrom] == cur_chr]  )

        df_part['cs'] = 0.0
        df_part['gen'] = ''
        df_part['tsd'] = 0
        df_part['updn'] = ''

        df_gtf_cur =  df_gtf.loc[df_gtf['seqname'] == gtf_chr].copy()

        df_gtf_cur = df_gtf_cur.reset_index(drop=True)

        df_gtf_cur_trans = df_gtf_cur.loc[df_gtf_cur['feature'] == 'transcript'].copy()

        idx_trans = 0

        N = 0

        wigFixchrFile = wig_dir + 'chr' + cur_chr + '.phastCons100way.wigFix'

        f = open(wigFixchrFile, 'r')
        last_pos = 0
        step = 0

        for index, row in df_part.iterrows():
            pos = row[columns.pos]

            #print(N)

            N = N + 1

            idx_trans, gen, dist, strand = ju.calc_gen_markup(pos, idx_trans, df_gtf_cur, df_gtf_cur_trans)
            valid, last_pos, step, score = ju.calc_cons_scores(pos, f, last_pos, step)

            df_part.at[index, 'cs'] = score
            df_part.at[index, 'gen'] = gen
            df_part.at[index, 'tsd'] = dist
            df_part.at[index, 'updn'] = strand

#        print(df_part.shape[0])
#        df_part = df_part.drop(df_part[df_part.cs < 0].index)
#        print(df_part.shape[0])

        if total_chr == 0:
            df_part.to_csv(output, sep='\t', index=False)
        else:
            df_part.to_csv(output, sep='\t', index=False, header=False, mode='a')

        gc.collect()

        total_chr = total_chr + 1

    print(datetime.datetime.now())
