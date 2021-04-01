import pandas as pd
import numpy as np
import configparser
import os
from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pywt


from gtfparse import read_gtf


import datetime
import gc

from argparse import ArgumentParser

class AlgParams:
    reverse_mode = False
    flexible_endings = False
    purine_orig_len = 3
    hydro_orig_len = 1
    purine_len = 1
    hydro_len = 1

alg_params = AlgParams()

pairs_p = pd.DataFrame( {'x' : pd.Series( [], dtype=int), 'y' : pd.Series( [], dtype=int)})
pairs_h = pd.DataFrame( {'x' : pd.Series( [], dtype=int), 'y' : pd.Series( [], dtype=int)})

pairs_p2 = pd.DataFrame( {'x' : pd.Series( [], dtype=int), 'y' : pd.Series( [], dtype=int)})
pairs_h2 = pd.DataFrame( {'x' : pd.Series( [], dtype=int), 'y' : pd.Series( [], dtype=int)})




def reset_df_data():
    global pairs_p
    global pairs_h
    global pairs_p2
    global pairs_h2

    pairs_p = pd.DataFrame({'x': pd.Series([], dtype=int), 'y': pd.Series([], dtype=int)})
    pairs_h = pd.DataFrame({'x': pd.Series([], dtype=int), 'y': pd.Series([], dtype=int)})
    pairs_p2 = pd.DataFrame({'x': pd.Series([], dtype=int), 'y': pd.Series([], dtype=int)})
    pairs_h2 = pd.DataFrame({'x': pd.Series([], dtype=int), 'y': pd.Series([], dtype=int)})



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
                if (i - start >= alg_params.purine_orig_len): # min length
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
                if (i - start >= alg_params.hydro_orig_len): # min length
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


def check_seq(seq, a_range, b_range, pairs, min_len = 0):
    seq_a = seq[a_range.start:a_range.end]
    seq_b = seq[b_range.start:b_range.end]

    if a_range.len >= min_len:
        if (seq_a == seq_b):
            pairs = set_matrix_value(a_range, b_range, pairs)
        elif( alg_params.reverse_mode and seq_a == seq_b[::-1] ):
            pairs = set_matrix_value_reverse(a_range, b_range, pairs)


    return pairs


def set_range(a_range, b_range, pairs, len = 0):
    if a_range.len >= len:
        pairs = set_matrix_value(a_range, b_range, pairs)
        pairs = set_matrix_value_reverse(a_range, b_range, pairs)

    return pairs

def get_helper_pts(len):
    x = []
    y_central = []

    y_above = []
    y_above2 = []
    y_above3 = []

    y_below = []
    y_below2 = []
    y_below3 = []

    for i in range(0, len):
        x.append(i)
        y_central.append(i)
        y_above.append(i + 75)
        y_above2.append(i + 150)
        y_above3.append(i + 300)

        y_below.append(i - 75)
        y_below2.append(i - 150)
        y_below3.append(i - 300)


    central_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_central, dtype=int)})
    above_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_above, dtype=int)})
    above2_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_above2, dtype=int)})
    above3_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_above3, dtype=int)})

    below_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_below, dtype=int)})
    below2_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_below2, dtype=int)})
    below3_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_below3, dtype=int)})

    return central_pts, above_pts, above2_pts, above3_pts, below_pts, below2_pts, below3_pts


def get_axes(len):
    x = []
    y = []

    x_zero = []
    y_zero = []


    for i in range(0, len):
        x.append(i)
        y.append(i)

        x_zero.append( int(len / 2) )
        y_zero.append( int(len / 2) )



    x_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y_zero, dtype=int)})
    y_pts = pd.DataFrame({'x': pd.Series(x_zero, dtype=int), 'y': pd.Series(y, dtype=int)})

    return x_pts, y_pts

def get_tandem_points(tandem_seq):
    x = []
    y = []
    for i in range(0, len(tandem_seq) ):
        if tandem_seq[i].islower():
            for k in range(0, len(tandem_seq), 8):
                x.append(i)
                y.append(k)
                x.append(k)
                y.append(i)

    tandem_pts = pd.DataFrame({'x': pd.Series(x, dtype=int), 'y': pd.Series(y, dtype=int)})

    return tandem_pts

def get_code_by_type(seq, type):

    if type == 'purine':
        anchor_seq = seq
        code1 = ['C', 'T']
        code2 = ['A', 'G']

    elif type == 'strong':
        anchor_seq = seq
        code1 = ['C', 'G']
        code2 = ['A', 'T']

    elif type == 'hydro':
        anchor_seq = to_hydro(seq)
        code1 = ['K', 'M']
        code2 = ['F', 'P']
    elif type == 'Gcode':
        anchor_seq = seq
        code1 = ['G']
        code2 = ['A', 'C', 'T']
    elif type == 'Acode':
        anchor_seq = seq
        code1 = ['A']
        code2 = ['G', 'C', 'T']
    else:
        assert(False)

    seq_np = np.zeros(len(anchor_seq), dtype=np.int32)

    for i in range(0, len(anchor_seq)):
        if anchor_seq[i] in code1:
            seq_np[i] = 1
        elif anchor_seq[i] in code2:
            seq_np[i] = -1
        elif anchor_seq[i] == 'N':
            seq_np[i] = 0
        else:
            assert (False)

    return seq_np


def make_wavelet(seq, type):
        seq_np = get_code_by_type(seq, type)
        idx_gaps = np.where(seq_np == 0)[0]
        seq_np_cur = seq_np.copy()

        seq_np_cur[idx_gaps] = np.random.binomial(n=1, p=0.5, size=len(idx_gaps))
        seq_np_cur[idx_gaps] = (seq_np_cur[idx_gaps] * 2) - 1

        scales = np.arange(1, 150)
        waveletname = 'morl'
        dt = 1

        [coefficients, frequencies] = pywt.cwt(seq_np_cur - np.mean(seq_np_cur), scales, waveletname, dt)

        idxs_del = np.where(frequencies > 0.5)[0] #according to Nyquist
        frequencies = np.delete(frequencies, idxs_del)
        coefficients = np.delete(coefficients, idxs_del, 0)

        power = (abs(coefficients)) ** 2
        period = [ round( 1 / x) for x in frequencies ]

        if(len(seq_np) > 100000):
            chrom_poses = [ round(x / 1000000, 1) for x in range(0, len(seq_np)) ]
            chrom_label = "chrom pos (Mb)"
        else:
            chrom_poses = [ int(x)  for x in range(0, len(seq_np)) ]
            chrom_label = "chrom pos (bp)"


        heatmap_pd = pd.DataFrame(
            data=power,  # values
            index=period,
            columns=chrom_poses)

        return heatmap_pd

def test():
    data_url = 'http://bit.ly/2cLzoxH'
    gapminder = pd.read_csv(data_url)
    print( np.unique( gapminder['continent'] ))

def test_visual():
    for k in range(0, 3):
        fig = plt.figure()
        ax = fig.add_subplot(111)  # The big subplot
        ax.set_aspect('equal', adjustable='box')

        plt.sca(ax)

        pairs_all = pd.DataFrame({'x': pd.Series(np.zeros(500), dtype=int), 'y': pd.Series(np.zeros(500), dtype=int)})
        sns.scatterplot(data=pairs_all, x="x", y="y", ax=ax, s=0.15)
        # plt.legend([], [], frameon=False)

        ax.legend(loc='upper left', markerscale=0.2, bbox_to_anchor=(1.04, 1), fontsize=2)
        ax.set_xlim(0, 100)
        ax.set_ylim(0, 100)

        # ax.set_xticklabels(ax.get_xticklabels(), fontsize = 4)
        # ax.set_yticklabels(ax.get_yticklabels(), fontsize = 4)

        plt.xticks(np.arange(0, 100 + 1, 5), fontsize=2 )
        plt.yticks(np.arange(0, 100 + 1, 5), fontsize=2 )

        ax.set_xlabel("")
        ax.set_ylabel("")

        fig.tight_layout()

        output_png = '{}_{}.png'.format('file1', k)
        plt.savefig(output_png, dpi=2400, bbox_inches='tight')
        fig.clear()
        plt.close(fig)

        fig = plt.figure()
        ax_w = fig.add_subplot(111)  # The big subplot
        # ax_w.set_aspect('equal', adjustable='box')

        plt.sca(ax_w)

        sns.set(font_scale=0.5)
        sns.heatmap(np.zeros((100, 100)), cbar_kws={'label': 'Energy'}, ax=ax_w, cmap='plasma')

        ax_w.set_xlim(0, 100)
        plt.xticks(fontsize=3)
        plt.yticks(fontsize=3)

        # ax_w.set_title(description, fontweight='bold', fontsize=4)

        fig.tight_layout()

        output_png = '{}_{}.png'.format('file2', k)
        plt.savefig(output_png, dpi=2400, bbox_inches='tight')
        fig.clear()
        plt.close(fig)

        matplotlib.rc_file_defaults()

if __name__ == '__main__':
    #test()
    #test_visual()
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

    alg_params.purine_orig_len = config['DEFAULT'].getint('purine_orig_len')
    alg_params.hydro_orig_len = config['DEFAULT'].getint('hydro_orig_len')
    alg_params.purine_len = config['DEFAULT'].getint('purine_len')
    alg_params.hydro_len = config['DEFAULT'].getint('hydro_len')

    gtf_file = 'hg38.ncbiRefSeq.gtf'
    df_gtf = read_gtf(gtf_file)

    all_chr = ['1']

    use_trans = False

    for cur_chr in all_chr:

        print(cur_chr)
        num_pages = 1

        if use_trans:
            gtf_chr = 'chr' + cur_chr
            df_gtf_cur = df_gtf.loc[df_gtf['seqname'] == gtf_chr].copy()
            df_gtf_cur_trans = df_gtf_cur.loc[ (df_gtf_cur['feature'] == 'transcript') & (df_gtf_cur['strand'] == '+')  ].copy()
            df_trans_to_plot = df_gtf_cur_trans.sample(n=30)
            trans_starts = [6424776, 6785460, 7771296, 9539972, 9943475, 10210570, 11054589, 12019833, 15236561, 15758792]
            trans_starts = df_trans_to_plot.loc[:, 'start'].tolist()
            num_pages = len(trans_starts)

        chr_filename = 'chr' + cur_chr + '.fa'

        full_filename = os.path.join(chrom_dir, chr_filename)
        full_filename = 'Test_sequence2.fa'
        #full_filename = 'test2_fp_100.fa'
        #output = 'test2_fp_100.png'

        if not os.path.exists(full_filename):
            continue

        fasta_seqs = SeqIO.parse(open(full_filename), 'fasta')

        for fasta in fasta_seqs:
            name, sequence = fasta.id, str(fasta.seq)
            break


        common_start = 30000000
        step = 5000



        for k in range(0, num_pages ):
            reset_df_data()

            if use_trans:
                trans_start = trans_starts[k]

                start_cur = int(trans_start - step / 2)
                end_cur = int(trans_start + step / 2)
                chunk_sequence = sequence[start_cur:end_cur]
            else:
                start_cur = 0
                end_cur = len(sequence)
                chunk_sequence = sequence

            print(start_cur, end_cur)
            #chunk_sequence = 'NnnnnnnnCaaaaaatCgatatatatCaaaaaatGgtgnnnnGtctctcaNnnnggggGtctctcg'
            #chunk_sequence = chunk_sequence.upper()

            #chunk_sequence = 'ATCTCAGTCCTGGGGAAGCTCATGATGGAACCTTGTCCAGGGGGTTCTAATCTCAGTCCTGGGGAAGCTCATGATGGAACCTTGTCCAGGGGGTTCTAATCTCAGTCCTGGGGAAGCTCATGATGGAACCTTGTCCAGGGGGTTCTAATCTCAGTCCTGGGGAAGCTCATGATGGAACCTTGTCCAGGGGGTTCTAATCTCAGTCCTGGGGAAGCTCATGATGGAACCTTGTCCAGGGGGTTCTA'
            tandem_sequence = chunk_sequence
            chunk_sequence = chunk_sequence.upper()

            purin_seq = to_purin(chunk_sequence)
            hydro_seq = to_hydro(chunk_sequence)

            purine_list = find_purine_chains(purin_seq)
            purine_list = sorted(purine_list)

            hydro_list = find_hydro_chains(hydro_seq)
            hydro_list = sorted(hydro_list)


            print(datetime.datetime.now())
            print('before calc')

            for cur_group in range (0, len(purine_list) - 1):
                for i in range(cur_group + 1, len(purine_list) ):
                    if purine_list[cur_group].len == purine_list[i].len:
                        pairs_p = check_seq(chunk_sequence, purine_list[cur_group], purine_list[i], pairs_p)
                        pairs_p2 = set_range(purine_list[cur_group], purine_list[i], pairs_p2, alg_params.purine_len)
                    else:
                        break


            for cur_group in range (0, len(hydro_list) - 1):
                for i in range(cur_group + 1, len(hydro_list) ):
                    if hydro_list[cur_group].len == hydro_list[i].len:
                        pairs_h = check_seq(chunk_sequence, hydro_list[cur_group], hydro_list[i], pairs_h)
                        pairs_h2 = check_seq(hydro_seq, hydro_list[cur_group], hydro_list[i], pairs_h2, alg_params.hydro_len)
                    else:
                        break

            print(datetime.datetime.now())
            print('after calc intervals')

            x_ax_pts, y_ax_pts = get_axes(len(chunk_sequence))
            central_pts, above_pts, above2_pts, above3_pts, below_pts, below2_pts,  below3_pts = get_helper_pts(len(chunk_sequence) )

            print(datetime.datetime.now())
            print('after calc lines')
            tandem_pts = get_tandem_points(tandem_sequence)

            print(datetime.datetime.now())
            print('after calc tandems')

            pairs_common = pd.merge(pairs_p, pairs_h, how="inner", on=["x", "y"])

            tandem_pts['type'] =   'Tandems'
            x_ax_pts['type'] = 'X_axis'
            y_ax_pts['type'] = 'Y_axis'

            central_pts['type'] = 'Diagonal'
            above_pts['type'] =   'Above 75'
            above2_pts['type'] =  'Above 150'
            above3_pts['type'] =  'Above 300'

            below_pts['type'] =   'Below 75'
            below2_pts['type'] =  'Below 150'
            below3_pts['type'] =  'Below 300'

            pairs_h2['type'] = 'Hydrocode by Hydrocode. ' + 'MinLen: {}'.format(alg_params.hydro_len)
            pairs_p2['type'] = 'Purine by Length. ' + 'MinLen: {}'.format(alg_params.purine_len)
            pairs_p['type'] = 'Purine by Nucleotides. ' + 'MinLen: {}'.format(alg_params.purine_orig_len)
            pairs_h['type'] = 'Hydrocode by Nucleotides. ' + 'MinLen: {}'.format(alg_params.hydro_orig_len)
            pairs_common['type'] = 'Common P & H by Nucleotides'

            pairs_all = pd.concat([tandem_pts, x_ax_pts, y_ax_pts, central_pts, above_pts, above2_pts, above3_pts, below_pts, below2_pts, below3_pts, \
                                   pairs_h2, pairs_p2, pairs_p, pairs_h, pairs_common])

            print(datetime.datetime.now())
            print('after concat. before plot')

            palette = []

            if tandem_pts.shape[0] > 0:
                palette.append('paleturquoise')

            palette.append('paleturquoise')
            palette.append('paleturquoise')

            palette.append('olivedrab')
            palette.append('lightgreen')
            palette.append('lightgreen')
            palette.append('lightgreen')
            palette.append('lightgreen')
            palette.append('lightgreen')
            palette.append('lightgreen')

            if pairs_h2.shape[0] > 0:
                palette.append('lightsteelblue')

            if pairs_p2.shape[0] > 0:
                palette.append('navajowhite')

            if pairs_p.shape[0] > 0:
                palette.append('brown')

            if pairs_h.shape[0] > 0:
                palette.append('royalblue')

            if pairs_common.shape[0] > 0:
                palette.append('black')

            fig = plt.figure()

            ax = fig.add_subplot(111)  # The big subplot
            ax.set_aspect('equal', adjustable='box')


            plt.sca(ax)

            sns.scatterplot(data=pairs_all, x="x", y="y", hue="type",  ax=ax, s=0.15,  palette=palette, legend='full')
            #plt.legend([], [], frameon=False)

            ax.legend(loc='upper left', markerscale=0.2, bbox_to_anchor=(1.04, 1), fontsize=2)
            ax.set_xlim(0, len(chunk_sequence))
            ax.set_ylim(0, len(chunk_sequence))

            #ax.set_xticklabels(ax.get_xticklabels(), fontsize = 4)
            #ax.set_yticklabels(ax.get_yticklabels(), fontsize = 4)

            plt.xticks(np.arange(0, len(chunk_sequence)+1, round( len(chunk_sequence) / 20) ), fontsize=2 )
            plt.yticks(np.arange(0, len(chunk_sequence)+1, round( len(chunk_sequence) / 20) ), fontsize=2)

            ax.set_xlabel("")
            ax.set_ylabel("")


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
            fig.clear()
            plt.close(fig)

            print(datetime.datetime.now())
            print('before wavelet')

            code_types = ['purine', 'strong', 'hydro', 'Gcode', 'Acode']

            for m in range(0, len(code_types) ):
                heatmap_wavelets = make_wavelet(chunk_sequence, code_types[m])


                fig = plt.figure()
                ax_w = fig.add_subplot(111)  # The big subplot
                #ax_w.set_aspect('equal', adjustable='box')

                plt.sca(ax_w)

                sns.set(font_scale=0.5)
                sns.heatmap(heatmap_wavelets, cbar_kws={'label': 'Energy'}, ax=ax_w, cmap = 'plasma')

                ax_w.set_xlim(0, len(chunk_sequence))
                plt.xticks(fontsize=3)
                plt.yticks(fontsize=3)

                #ax_w.set_title(description, fontweight='bold', fontsize=4)

                fig.tight_layout()

                output_png = '{}_{}_wavelet_type_{}.png'.format(output, k, code_types[m])
                plt.savefig(output_png, dpi=2400, bbox_inches='tight')
                fig.clear()
                plt.close(fig)

                matplotlib.rc_file_defaults()


    print(datetime.datetime.now())


