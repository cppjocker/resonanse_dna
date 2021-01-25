from gtfparse import read_gtf
import pandas as pd
import numpy as np
import configparser

import datetime

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
    else:
        pass

def calc_cons_scores(pos_SNP, f, last_pos, step):
    score = -1.0

    while last_pos < pos_SNP:

        line = f.readline()

        if not line:
            return False, last_pos, step, -1

        last_pos = last_pos + step

        parts = line.split()
        assert( len(parts) == 4 or len(parts) == 1 )

        if len(parts) == 4:
            assert( parts[0] == 'fixedStep')

            chrom = parts[1]
            start = parts[2]
            new_step = parts[3]

            start_spl = start.split('=')
            new_step_spl = new_step.split('=')

            assert(chrom.split('=')[0] == 'chrom')
            assert(start_spl[0] == 'start')
            assert(new_step_spl[0] == 'step')

            last_pos = int( start_spl[1] )
            step = int(new_step_spl[1])

            line = f.readline()
            assert (line)

            parts = line.split()
            assert (len(parts) == 1)

        assert(len(parts) == 1)

        score = float(line)

    if last_pos > pos_SNP:
        return False, last_pos, step, -1

    return True, last_pos, step, score





def calc_gen_markup(pos, idx_trans, df_gtf_cur, df_gtf_cur_trans):
    found_intersec = False
    return_pos = 0

    has_exon = False
    has_intron = 0
    has_3utr = 0
    has_5utr = 0

    has_itrgen = False

    idx_trans_start = idx_trans

    while True:
        if idx_trans >= df_gtf_cur_trans.shape[0] - 1:
            break

        start_trans = df_gtf_cur_trans.iat[idx_trans, 3]
        end_trans = df_gtf_cur_trans.iat[idx_trans, 4]

        next_trans_start = df_gtf_cur_trans.iat[idx_trans + 1, 3]

        assert (start_trans <= next_trans_start)

        if start_trans > pos:
            break

        if (not found_intersec) and (start_trans <= pos <= end_trans):
            found_intersec = True
            return_pos = idx_trans

        if start_trans <= pos <= end_trans:
            orig_index_start = df_gtf_cur_trans.index[idx_trans]
            orig_index_end = df_gtf_cur_trans.index[idx_trans + 1]

            has_exon_cur = False
            has_3utr_cur = 0
            has_5utr_cur = 0

            for i in range(orig_index_start, orig_index_end):
                cur_feature = df_gtf_cur.iat[i, 2]
                cur_feature_start = df_gtf_cur.iat[i, 3]
                cur_feature_end = df_gtf_cur.iat[i, 4]

                if not (cur_feature in {'CDS', '3UTR', '5UTR'}):
                    continue

                if not (cur_feature_start <= pos <= cur_feature_end):
                    continue

                if cur_feature == 'CDS':
                    has_exon_cur = True

                if cur_feature == '3UTR':
                    has_3utr_cur = 1

                if cur_feature == '5UTR':
                    has_5utr_cur = 1

            has_exon = has_exon or has_exon_cur

            if has_3utr_cur == 1:
                has_3utr = has_3utr_cur

            if has_5utr_cur == 1:
                has_5utr = has_5utr_cur

            if (not has_exon_cur) and (has_3utr_cur + has_5utr_cur == 0):
                has_intron = True

        idx_trans = idx_trans + 1

    idx_trans_end = idx_trans

    idx_trans_start = max(idx_trans_start - 20, 0)
    idx_trans_end = min(idx_trans_end + 20, df_gtf_cur_trans.shape[0])

    dist = 1000000000000
    strand = 'xx'

    for i in range(idx_trans_start, idx_trans_end):
        if df_gtf_cur_trans.iat[i, 6] == '+':  # strand
            cur_dist = abs(pos - df_gtf_cur_trans.iat[i, 3])  # start
            if pos < df_gtf_cur_trans.iat[i, 3]:
                cur_strand = 'up'
            else:
                cur_strand = 'down'
        elif df_gtf_cur_trans.iat[i, 6] == '-':
            cur_dist = abs(pos - df_gtf_cur_trans.iat[i, 4])  # end
            if pos > df_gtf_cur_trans.iat[i, 3]:
                cur_strand = 'up'
            else:
                cur_strand = 'down'
        else:
            assert (False)

        if cur_dist < dist:
            dist = cur_dist
            strand = cur_strand

    if found_intersec:
        idx_trans = return_pos

    weird_sum = has_intron + has_3utr + has_5utr

    if (not has_exon) and (weird_sum == 0):
        has_itrgen = True

    gen = 'unk'

    if has_itrgen:
        gen = 'igc'
    elif has_exon:
        gen = 'CDS'
    elif has_intron and weird_sum == 1:
        gen = 'tro'
    elif has_3utr and weird_sum == 1:
        gen = 'tpu'
    elif has_5utr and weird_sum == 1:
        gen = 'fpu'

    return idx_trans, gen, dist, strand

if __name__ == '__main__':
    print(datetime.datetime.now())

    df_gtf = read_gtf("hg19.ncbiRefSeq.gtf")

    config = configparser.ConfigParser()

    config.read('config_gen.ini')

    chrom_dir = config['DEFAULT']['chrom_dir']
    wig_dir = config['DEFAULT']['wig_dir']

    snp_table_file = config['DEFAULT']['input_file']
    output = config['DEFAULT']['output_file']

    pd_header = pd.read_csv(snp_table_file, sep='\t', header=0, nrows=3)
    determineColumns(pd_header)

    df_snp = pd.read_csv(snp_table_file, sep='\t', header=0, dtype={columns.chrom: str})

    df_snp = df_snp.dropna(subset=[columns.pos])

    if (df_snp[columns.chrom].dtype == np.int64) or (df_snp[columns.chrom].dtype == np.int32):
        df_snp[columns.chrom].dtype

    df_snp = df_snp.astype({columns.chrom: str, columns.pos: int})

    print(set(df_snp[columns.chrom]))


    all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
               '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    #all_chr = ['18']

    df_snp['cs'] = 0.0
    df_snp['gen'] = ''
    df_snp['tsd'] = 0
    df_snp['updn'] = ''

    df_write = pd.DataFrame(columns=df_snp.columns)


    for cur_chr in all_chr:
        print(cur_chr)

        gtf_chr = 'chr' + cur_chr

        df_part = df_snp.loc[df_snp[columns.chrom] == cur_chr].copy()

        df_gtf_cur =  df_gtf.loc[df_gtf['seqname'] == gtf_chr].copy()

        df_gtf_cur = df_gtf_cur.reset_index(drop=True)

        df_gtf_cur_trans = df_gtf_cur.loc[df_gtf_cur['feature'] == 'transcript'].copy()

        idx_trans = 0

        df_part_write = df_part.copy()

        N = 0

        wigFixchrFile = wig_dir + 'chr' + cur_chr + '.phastCons100way.wigFix'

        f = open(wigFixchrFile, 'r')
        last_pos = 0
        step = 0

        for index, row in df_part.iterrows():
            pos = row[columns.pos]

            #print(N)

            N = N + 1


            idx_trans, gen, dist, strand = calc_gen_markup(pos, idx_trans, df_gtf_cur, df_gtf_cur_trans)
            valid, last_pos, step, score = calc_cons_scores(pos, f, last_pos, step)

            df_part_write.at[index, 'cs'] = score
            df_part_write.at[index, 'gen'] = gen
            df_part_write.at[index, 'tsd'] = dist
            df_part_write.at[index, 'updn'] = strand

        df_write = pd.concat([df_write, df_part_write])

    print(df_write.shape[0])
    df_write_drop = df_write.drop(df_write[df_write.cs < 0].index)
    print(df_write_drop.shape[0])

    #df_part_write = df_part_write.iloc[0:1000]

    df_write.to_csv(output, sep='\t', index=False)

    print(datetime.datetime.now())
