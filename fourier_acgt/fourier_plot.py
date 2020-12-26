import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import math

import configparser
from argparse import ArgumentParser
import argparse
import os

def dir_path(path):
    if os.path.isdir(path):
        return path
    else:
        raise NotADirectoryError

if __name__ == '__main__':
    config = configparser.ConfigParser()

    parser = ArgumentParser()
    parser.add_argument("-conf", "--config", help="config ini file", type=argparse.FileType('r'), required=False)

    parser.add_argument("-in", "--input",
                        help="heatmap file", type=argparse.FileType('r'))

    parser.add_argument("-out", "--output",
                        help="png file", type=argparse.FileType('w'))

    parser.add_argument("-indir", "--input_dir",
                        help="input directory of fasta file", type = dir_path)

    parser.add_argument("-outdir", "--output_dir",
                        help="input directory of fasta file")

    parser.add_argument("-perc", "--percentile", help = "percentile of visualization")

    args = parser.parse_args()

    if args.config is not None:
        config.read(args.conf)
    else:
        config.read('config_v.ini')

    if args.input is not None:
        config['DEFAULT']['heatmap_file'] = args.input

    if args.output is not None:
        config['DEFAULT']['png_file'] = args.out

    if args.input_dir is not None:
        config['DEFAULT']['input_dir'] = args.input_dir

    if args.output_dir is not None:
        config['DEFAULT']['output_dir'] = args.output_dir

    if args.percentile is not None:
        config['VISUAL'] ['percentile'] = args.percentile

    if not os.path.exists(config['DEFAULT']['output_dir'] ) :
        os.mkdir(config['DEFAULT']['output_dir'])



    #"heatmap.csv"

    heatmap_file = os.path.join(config['DEFAULT']['input_dir'],  config['DEFAULT']['heatmap_file'])
    png_file = os.path.join(config['DEFAULT']['output_dir'],  config['DEFAULT']['png_file'])

    percentile = config['VISUAL'].getint('percentile')
    n_chr_ranges = config['VISUAL'].getint('chr_chunks')

    heatmap_mat = pd.read_csv(heatmap_file, sep=";", skiprows=1, index_col=0, header=0)

    with open(heatmap_file, 'r') as fin:
        description = fin.readline()

    print("Loaded")

    f_ranges = [(0.0, 0.002), (0.002, 0.004), (0.004, 0.007), (0.007, 0.01), (0.01, 0.02), (0.02, 0.03), (0.03, 0.04), (0.04, 0.05), (0.05, 0.1), (0.1, 0.2), (0.2, 0.3), (0.3, 0.4), (0.4, 0.5)  ]

    total_chr_chunks = heatmap_mat.shape[1]

    visual_chunk = math.floor( total_chr_chunks / n_chr_ranges )

    fig, axs = plt.subplots( len(f_ranges), n_chr_ranges, sharex='col', sharey='row')

    plt.suptitle(description + 'Percentile : {}.'.format(percentile), fontsize = 7)
   # plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)


    for k in range(0, len(f_ranges ) ):
        for j in range(0, n_chr_ranges):
            chr_start = j * visual_chunk
            chr_end = j * visual_chunk + visual_chunk

            f_range = f_ranges[k]
            f_start = f_range[0]
            f_end = f_range[1]


            heatmap_mat_chunk = heatmap_mat.loc[(heatmap_mat.index > f_start) & (heatmap_mat.index <= f_end)].iloc[:, chr_start:chr_end ].copy()

            columns_num = [round(float(x) / 1000000, 1 ) for x in heatmap_mat_chunk.columns]

            #change matrix by percentile logic

            for i_hmap in range(0, heatmap_mat_chunk.shape[0]):
                perc_xx = np.percentile(heatmap_mat_chunk.iloc[i_hmap, :], percentile)
                idxs_more_perc_xx = np.where(heatmap_mat_chunk.iloc[i_hmap, :] > perc_xx)[0]
                heatmap_mat_chunk.iloc[i_hmap, idxs_more_perc_xx] = perc_xx
                pass

            print(heatmap_mat_chunk.shape)

            font_scale = 0.2

    #            if(heatmap_mat_chunk.shape[0] < 40):
    #                font_scale = 0.2

            round_val = 1


            indexes = [round(1 / x, round_val ) for x in heatmap_mat_chunk.index]

            min_y = min( indexes)
            max_y = max(indexes)

            yticks = int( heatmap_mat_chunk.shape[0] / 5)


            if (max_y - min_y) > 5:
                round_val = 0

            print(yticks, round_val)

            if (round_val == 0):
              indexes = [ int(x)  for x in indexes]

            heatmap_mat_chunk.index = indexes
            heatmap_mat_chunk.columns = columns_num


            sns.set(font_scale=font_scale)
            if axs.ndim == 2:
                sns_plt = sns.heatmap(heatmap_mat_chunk, cbar_kws={'label': 'Energy'}, cmap = 'plasma', yticklabels=yticks, ax = axs[k, j])
            else:
                sns_plt = sns.heatmap(heatmap_mat_chunk, cbar_kws={'label': 'Energy'}, cmap = 'plasma', yticklabels=yticks, ax = axs[k])

            #sns.color_palette("turbo", as_cmap=True)

            #sns_plt.set_yticks(  np.linspace(min_y, max_y,  5 ) )


            sns_plt.set_xticklabels( sns_plt.get_xmajorticklabels(), fontsize = 3 )
            sns_plt.set_yticklabels( sns_plt.get_ymajorticklabels(), fontsize = 3)


            sns_plt.set_xlabel("chrom pos (Mb)", fontsize=3)
            sns_plt.set_ylabel("period (bp)", verticalalignment='center', fontsize=3)


    plt.savefig(png_file, dpi = 1200 )

    plt.close()
