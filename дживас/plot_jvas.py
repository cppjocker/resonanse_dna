import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import os
from Bio import SeqIO






def case_AT_get():
    code = 'AT'

    fields = ['other_allele', 'effect_allele', 'effect_allele_frequency']

    sh_left = code + '_shoulder_l'
    sh_right = code + '_shoulder_r'
    fields.append(sh_left)
    fields.append(sh_right)
    fields.append('chromosome')
    fields.append('variant_id')
    fields.append('base_pair_location')

    df = pd.read_csv('Howard_hydro_hd07_gen.tsv', usecols=fields, sep='\t', header=0)
    df = df.reindex(columns=fields)

    df = df.astype({'chromosome': str, 'base_pair_location': int})

    df = df.loc[df.iloc[:, 3] >= 0, :]
    df_filter = df.loc[(df.iloc[:, 0] == 'G') | (df.iloc[:, 1] == 'G'), :] #AT->GT only
    df_filter['sum_shoulders'] = (df_filter[sh_left] + df_filter[sh_right])

    df_filter.sort_values(by=['sum_shoulders'], inplace=True)

    n_bins = 20
    amount_in_bin = df_filter.shape[0] // n_bins
    start_bin = (19 * amount_in_bin)
    end_bin = (19 * amount_in_bin + amount_in_bin)

    df_filter_end = df_filter.iloc[start_bin:end_bin, :]
    df_filter_end['seq'] = ''

    all_chr = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11',
               '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y']

    #all_chr = '1'
    idx_chr = 0
    for cur_chr in all_chr:
        print(cur_chr)

        df_cur_chr = df_filter_end.loc[ df_filter_end['chromosome'] == cur_chr]
        print(df_cur_chr.shape[0])

        chr_filename = 'chr' + cur_chr + '.fa'
        chr_full_filename = os.path.join("../hg19", chr_filename)

        print(chr_full_filename)

        fasta_seqs = SeqIO.parse(open(chr_full_filename), 'fasta')

        for fasta in fasta_seqs:
            name, sequence = fasta.id, str(fasta.seq)
            break

        sequence = sequence.lower()
        for i in df_cur_chr.index:
            pos = df_cur_chr.loc[i, 'base_pair_location']
            pos = pos - 1

            left_shoulder_len = df_cur_chr.loc[i, sh_left] + 2
            right_shoulder_len = df_cur_chr.loc[i, sh_right] + 2

            print(pos, left_shoulder_len, right_shoulder_len)

            cur_seq = sequence[ (pos - left_shoulder_len) :  (pos + right_shoulder_len + 2) ]

            print(df_cur_chr.loc[i, 'chromosome'])

            cur_seq = list(cur_seq)

            cur_seq[0] = cur_seq[0].upper()
            cur_seq[1] = cur_seq[1].upper()

            cur_seq[-1] = cur_seq[-1].upper()
            cur_seq[-2] = cur_seq[-2].upper()

            cur_seq[left_shoulder_len] = cur_seq[left_shoulder_len].upper()
            cur_seq[left_shoulder_len + 1] = cur_seq[left_shoulder_len + 1].upper()

            cur_seq = "".join(cur_seq)

            print(cur_seq)

            df_cur_chr.loc[i, 'seq'] = cur_seq

        df_cur_chr = df_cur_chr.loc[:, ['sum_shoulders', 'seq']]

        if idx_chr == 0:
            df_cur_chr.to_csv('Howard_GT_example.tsv', sep='\t', index=False)
        else:
            df_cur_chr.to_csv('Howard_GT_example.tsv', sep='\t', header = False, index=False, mode='a')

        idx_chr += 1


def case_purine_shoulder():
    shoulder_codes = ['AC', 'AT', 'GC', 'GT', 'CA', 'CG', 'TA', 'TG']
    shoulder_codes = ['AC', 'AT', 'GC', 'GT']

    #fields = ['other_allele', 'effect_allele', 'effect_allele_frequency',  'arm_AT', 'arm_AT_l', 'arm_AT_r'] #Howard
    #fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira

    #fields = ['A1', 'A2', 'freq2'] #freq2, GreatTit

    n_bins = 16

    out_df = pd.DataFrame( {'number_bin' : np.arange(n_bins) })
    filter_gen = 'igc_and_tro' #CDS, igc, tro
    start_cs = 0.0
    end_cs = 0.4
    cs_str = 'low'

    for code in shoulder_codes:
        fields = ['other_allele', 'effect_allele', 'effect_allele_frequency', 'r1_len', 'r2_len', 'y1_len', 'y2_len']  # Karlsson

        fields.append('gen')
        fields.append('cs')


        df = pd.read_csv('Howard_hydro_hd07_gen_fixcodes.tsv', usecols=fields, sep='\t', header=0)
        df = df.reindex(columns=fields)
        df = df.loc[ df.iloc[:, 3] >= 0, :]
        df = df.loc[df.loc[:, 'gen'] != 'CDS']
        df = df.loc[ (df.loc[:, 'cs'] >= start_cs) & (df.loc[:, 'cs'] <= end_cs)]

        #df = df.astype({code + '_AT_in': bool})

        main_letter = code[0]
        filters = [ 'A', 'C', 'G', 'T' ]
        filters.remove(main_letter)



        df_filter = df.loc[ ( (df.iloc[:, 0] == code[0]) & (df.iloc[:, 1] ==  code[1]) ) |
                            ( (df.iloc[:, 1] == code[0]) & (df.iloc[:, 0] ==  code[1])), :]

        def len_purine(row):
            if code[0] in ['A', 'G']:
                if row['other_allele'] in ['A', 'G']:
                    return row['r1_len']
                else:
                    return row['r2_len']
            else:
                if row['other_allele'] in ['C', 'T']:
                    return row['y1_len']
                else:
                    return row['y2_len']

        def freq_purine(row):
            if code[0] in ['A', 'G']:
                if row['other_allele'] in ['A', 'G']:
                    return 1 - row['effect_allele_frequency']
                else:
                    return row['effect_allele_frequency']
            else:
                if row['other_allele'] in ['C', 'T']:
                    return 1 - row['effect_allele_frequency']
                else:
                    return row['effect_allele_frequency']


        df_filter['purine_len'] = df_filter.apply(lambda row : len_purine(row), axis = 1)
        df_filter['purine_freq'] = df_filter.apply(lambda row : freq_purine(row), axis = 1)

        df_filter.sort_values(by=['purine_len'], inplace=True)

        bins = []
        freqs = []
        amount_in_bin = df_filter.shape[0] // n_bins

       # print(df_filter['purine_len'])

        for i in range(0, n_bins):
            start_bin = (i * amount_in_bin)
            end_bin = (i * amount_in_bin + amount_in_bin)


            #bin_center = np.median(df_filter.iloc[start_bin:end_bin, 8 ])
            bin_center = i
            #bin_center = int( start_bin+ amount_in_bin / 2 )
            #bin_center = np.median( df_filter.iloc[ start_bin:end_bin, -2 ].values  )
            #bin_center = np.median( 0.5 * ( df_filter.iloc[start_bin:end_bin, 3 ].values + df_filter.iloc[start_bin:end_bin, 4 ].values ) )

            cur_freqs = df_filter.iloc[start_bin:end_bin, -1]

            bins.append(bin_center)
            freqs.append(np.mean(cur_freqs) )

        names = '{0}>{1}'.format(code[0], code[1])

        plt.plot(bins, freqs, '-o', label = names)

        cur_df = pd.DataFrame({'bin_SNP_{}'.format(code) : bins, 'freq_SNP_{}'.format(code): freqs})
        out_df = pd.concat( [out_df, cur_df], axis=1)

    #plt.axvline(x=80000, c='r')
    plt.xlabel('Shoulder bin')
    plt.ylabel('allele frequency')
    plt.title('Purine Evolution.  Gen {}. CS {:.2f}-{:.2f} '.format(filter_gen, start_cs, end_cs))
    plt.legend(loc='best')
    plt.savefig('Howard_bins_purine_shoulder_evolution_average_gen_{}_cs_{}.png'.format( filter_gen, cs_str), dpi=600, bbox_inches='tight')
    plt.close()

    out_df.to_csv('Howard_bins_purine_shoulder_CS_{}_gen_{}.tsv'.format(cs_str, filter_gen), sep='\t', index=False)




def caseAT():
    shoulder_codes = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'GA', 'GC', 'TA']

    #fields = ['other_allele', 'effect_allele', 'effect_allele_frequency',  'arm_AT', 'arm_AT_l', 'arm_AT_r'] #Howard
    #fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira

    #fields = ['A1', 'A2', 'freq2'] #freq2, GreatTit

    n_bins = 6

    out_df = pd.DataFrame( {'number_bin' : np.arange(n_bins) })
    filter_gen = 'igc_and_tro' #CDS, igc, tro
    start_cs = 0.0
    end_cs = 0.4
    cs_str = 'low'

    for code in shoulder_codes:
        fields = ['other_allele', 'effect_allele', 'effect_allele_frequency']  # Karlsson

        sh_left = code + '_shoulder_l'
        sh_right = code + '_shoulder_r'
        fields.append(sh_left)
        fields.append(sh_right)
        fields.append('gen')
        fields.append('cs')
        fields.append(code + '_AT_in')


        df = pd.read_csv('Howard_hydro_hd07_gen_fixcodes.tsv', usecols=fields, sep='\t', header=0)
        df = df.reindex(columns=fields)
        df = df.loc[ df.iloc[:, 3] >= 0, :]
        df = df.loc[df.loc[:, 'gen'] != 'CDS']
        df = df.loc[ (df.loc[:, 'cs'] >= start_cs) & (df.loc[:, 'cs'] <= end_cs)]

        #df = df.astype({code + '_AT_in': bool})

        main_letter = code[0]
        filters = [ 'A', 'C', 'G', 'T' ]
        filters.remove(main_letter)

        overall_bins = []
        overall_freqs = []
        overall_names = []

        filters_AT_inside = ['b', 'u']

        for f_AT, f_AT_bool in zip(filters_AT_inside, [True, False] ):
            for filter in filters:

                df_filter = df.loc[ (df.iloc[:, 0] == filter) | (df.iloc[:, 1] == filter), :]
                df_filter = df_filter.loc[ (df_filter.iloc[:, -1] == f_AT_bool), :]


                df_filter['arm_AT'] = 0.5 * ( df_filter[sh_left] + df_filter[sh_right] )
                df_filter.sort_values(by=['arm_AT'], inplace=True)

                bins = []
                freqs = []
                amount_in_bin = df_filter.shape[0] // n_bins

                # all_freqs = [df_filter.iloc[x, 2] if df_filter.iloc[x, 1] == main_letter else 1 - df_filter.iloc[x, 2] for x in
                #              range(0, df_filter.shape[0])]
                #
                # fig = plt.figure()
                # ax = Axes3D(fig)
                # surf = ax.plot_trisurf(df_filter.iloc[:, 3], df_filter.iloc[:, 4], all_freqs, cmap=cm.jet, linewidth=0.1)
                # ax.set_xlabel('l1')
                # ax.set_ylabel('l2')
                # ax.set_zlabel('freq ' + main_letter)
                #
                # fig.colorbar(surf, shrink=0.5, aspect=5)
                # ax.set_title('Average. Denucl: {0} SNP: {1}->{2}'.format(code, main_letter, filter))
                # plt.savefig('Karlsson_3D_SNP_{}_{}_{}_Multic_AT2.png'.format(code, main_letter, filter), dpi=300)
                #
                # plt.close()
                # #plt.pause(0)


                for i in range(0, n_bins):
                    start_bin = (i * amount_in_bin)
                    end_bin = (i * amount_in_bin + amount_in_bin)

                    bin_center = np.median(df_filter.iloc[start_bin:end_bin, 8 ])
                    bin_center = i
                    #bin_center = int( start_bin+ amount_in_bin / 2 )

                    #bin_center = np.median( np.sqrt( df_filter.iloc[start_bin:end_bin, 3 ].values ) )

                    #bin_center = np.median( 0.5 * ( df_filter.iloc[start_bin:end_bin, 3 ].values + df_filter.iloc[start_bin:end_bin, 4 ].values ) )


                    cur_freqs = [ df_filter.iloc[x, 2]  if df_filter.iloc[x, 1] == main_letter else 1 - df_filter.iloc[x, 2] for x in range(start_bin, end_bin) ]


                    bins.append(bin_center)
                    freqs.append(np.mean(cur_freqs) )



                # plt.plot(bins, freqs, '-o')
                # plt.xlabel('Average shoulder')
                # plt.ylabel('allele AT frequency')
                # plt.title('Multic AT_2. SNP: {0}->{1}'.format(main_letter, filter))
                # plt.savefig('Great_tit_arm_evolution_SNP_{}_Average.png'.format(main_letter + filter), dpi=600, bbox_inches='tight')
                # plt.close()
                #
                # cur_df = pd.DataFrame({'bins' : bins, 'freqs' : freqs})
                #
                # cur_df.to_csv('evolution_arm_Great_tit_{}_Average.tsv'.format(main_letter + filter), sep='\t', index=False)

                overall_bins.append(bins)
                overall_freqs.append(freqs)
                overall_names.append('{0}>{1} {2}'.format(code, filter + code[1] , f_AT))

        for i in range( len(overall_freqs) ):
            bins = overall_bins[i]
            freqs =  overall_freqs[i]
            names = overall_names[i]

            plt.plot(bins, freqs, '-o', label = names)

            cur_df = pd.DataFrame({'bin_dinucl_{}_SNP_{}'.format(code, names) : bins, 'freq_dinucl_{}_SNP_{}'.format(code, names): freqs})

            out_df = pd.concat( [out_df, cur_df], axis=1)

        #plt.axvline(x=80000, c='r')
        plt.xlabel('Shoulder bin length number')
        plt.ylabel('allele {} frequency'.format(code))
        plt.title('Dinucl {}. Gen {}. CS {:.2f}-{:.2f} '.format(code, filter_gen, start_cs, end_cs))
        plt.legend(loc='best')
        plt.savefig('Howard_num_arm_evolution_all_SNP_Average_Dinucl_{}_gen_{}_cs_{}.png'.format(code, filter_gen, cs_str), dpi=600, bbox_inches='tight')
        plt.close()

    out_df.to_csv('Howard_num_all_SNP_all_dinucl_{}_CS_{}.tsv'.format(filter_gen, cs_str), sep='\t', index=False)





def case1():
    #fields = ['other_allele', 'effect_allele', 'effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    # fields = ['EAF_A1',  'm1', 'm2', 'p1', 'p2'] #Karlsson
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira
    fields = ['A1', 'A2', 'freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit

    df = pd.read_csv('Great_tit_hydro.tsv', usecols=fields, sep='\t', header=0)

    filters = [['A', 'G'], ['A','C'], ['A', 'T'], ['C','G'], ['C','T'], ['G', 'T']]

    #filters = [['A', 'G', 'C', 'T']]
    print(df.head)

    points = np.arange(0.5, 1.01, 0.05)
    bins = [0.5 * (points[i] + points[i + 1]) for i in range(0, len(points) - 1)]
    df_out = pd.DataFrame( {'bins' : bins} )


    for min_diff in [0, 2, 4, 6]:
        for filter in filters:

            counts_h = np.zeros( len(bins) )
            counts_p = np.zeros( len(bins) )

            sum_lens_p = np.zeros( len(bins) )
            sum_lens_h = np.zeros( len(bins) )


            for k in range(0, df.shape[0]):
                col_shift = 0
                if len(fields ) > 5:
                    a1 = df.iloc[k, 0]
                    a2 = df.iloc[k, 1]
                    col_shift = 2

                freq2 = df.iloc[k, 0 + col_shift]
                freq1 = 1 - freq2

                m1 = df.iloc[k, 1 + col_shift]
                m2 = df.iloc[k, 2 + col_shift]

                p1 = df.iloc[k, 3 + col_shift]
                p2 = df.iloc[k, 4 + col_shift]


                if (a1 not in filter) or (a2 not in filter):
                    continue

                freq = max(freq1, freq2)

                if  freq1 > freq2:
                    freq = freq1
                    diff_m = m1 - m2
                    diff_h = p1 - p2
                else:
                    freq = freq2
                    diff_m = m2 - m1
                    diff_h = p2 - p1

                for i in range(0, len(bins) ):
                    if( points[i] <= freq  and freq < points[i + 1] ):
                        if abs(diff_m) >= min_diff:
                            counts_p[i] += 1
                            sum_lens_p[i] += diff_m
                        if abs(diff_h) >= min_diff:
                            counts_h[i] += 1
                            sum_lens_h[i] += diff_h

            # print(bins)
            # print(counts_p)
            # print(counts_h)
            #
            # print(sum_lens_p)
            # print(sum_lens_h)


            sum_lens_p = [ sum_lens_p[i] / counts_p[i] for i in range(0, len(bins))]
            sum_lens_h = [ sum_lens_h[i] / counts_h[i] for i in range(0, len(bins))]


            plt.plot(bins, sum_lens_p)
            plt.xlabel('frquency bin')
            plt.ylabel('mean diff len')
            plt.title('purine. Min Diff: {} bp'.format(min_diff))
            plt.savefig('GreatTit_purine_evolution_diff_{}_filt_{}.png'.format(min_diff, filter[0] + filter[1]), dpi=600, bbox_inches='tight')
            plt.close()

            plt.plot(bins, sum_lens_h)
            plt.xlabel('frquency bin')
            plt.ylabel('mean diff len')
            plt.title('hydro. Min Diff: {} bp'.format(min_diff))
            plt.savefig('GreatTit_hydro_evolution_diff_{}_filt_{}.png'.format(min_diff, filter[0] + filter[1]), dpi=600, bbox_inches='tight')
            plt.close()

            cur_df = pd.DataFrame({'{}-{}_purine'.format(min_diff, filter[0] + filter[1]) : sum_lens_p, '{}-{}_hydro'.format(min_diff, filter[0] + filter[1]) : sum_lens_h})
            df_out = pd.concat([df_out, cur_df], axis=1)

    df_out.to_csv('effects_great_tit.tsv', sep='\t', index=False)


def case2():
    #fields = ['variant_id', 'effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard
    #fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard

    #fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang, LeeJJ, Mahajan
    #fields = ['EAF_A1',  'm1', 'm2', 'p1', 'p2'] #Karlsson
    #fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    #fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira

    fields = ['freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit
    #fields = ['CHROM', 'POS', 'freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit


    #fields = ['ExAC_MAF',  'm1', 'm2', 'p1', 'p2'] #Turcot


    #df = pd.read_csv('Mahajan_hydro.tsv', usecols=fields,sep='\t', header=0)
    df = pd.read_csv('Great_tit_hydro.tsv', usecols=fields, sep='\t', header=0)



    # df = df.loc[ (df['variant_id'] == 'rs4951862' ) |  (df['variant_id'] == 'rs4500250') ]
    # assert(df.shape[0] == 2)
    # df = df.iloc[:, 1:]

    #'Great_tit_hydro.tsv'
    #'Howard_test_hydro_ext.tsv'
    #'salmon_hydro_test.tsv'

    print(df.head)

    points = np.arange(0.5, 1.01, 0.05)
    bins = [0.5 * (points[i] + points[i + 1]) for i in range(0, len(points) - 1)]
    df_out = pd.DataFrame( {'bins' : bins} )


    for range_len in [ [1, 1], [2, 2], [3, 3], [4, 5], [6, 10], [11, 15], [16, 10000] ]:
    #for range_len in [[1, 10000]]:

            range_s = range_len[0]
            range_e = range_len[1]


            counts_h = np.zeros( len(bins) )
            counts_p = np.zeros( len(bins) )

            counts_p_g = np.zeros( len(bins) )
            counts_p_e = np.zeros( len(bins) )
            counts_p_l = np.zeros( len(bins) )

            counts_h_g = np.zeros( len(bins) )
            counts_h_e = np.zeros( len(bins) )
            counts_h_l = np.zeros( len(bins) )

            sum_lens_p = np.zeros( len(bins) )
            sum_lens_h = np.zeros( len(bins) )


            for k in range(0, df.shape[0]):
                col_shift = 0

                freq2 = df.iloc[k, 0]
                freq1 = 1 - freq2

                m1 = df.iloc[k, 1]
                m2 = df.iloc[k, 2]
                p1 = df.iloc[k, 3]
                p2 = df.iloc[k, 4]

                min_m = min(m1, m2)
                min_p = min(p1, p2)


                freq = max(freq1, freq2)

                if  freq1 > freq2:
                    freq = freq1
                    diff_m = m1 - m2
                    diff_h = p1 - p2
                else:
                    freq = freq2
                    diff_m = m2 - m1
                    diff_h = p2 - p1

                for i in range(0, len(bins) ):
                    if( points[i] <= freq  and freq < points[i + 1] ):
                        if range_s <= min_m and min_m <= range_e:
                            counts_p[i] += 1
                            sum_lens_p[i] += diff_m

                            if diff_m > 0:
                                counts_p_g[i] += 1
                            elif diff_m == 0:
                                counts_p_e[i] += 1
                            else:
                                counts_p_l[i] += 1

                        if range_s <= min_p and min_p <= range_e:
                            counts_h[i] += 1
                            sum_lens_h[i] += diff_h

                            if diff_h > 0:
                                counts_h_g[i] += 1
                            elif diff_h == 0:
                                counts_h_e[i] += 1
                            else:
                                counts_h_l[i] += 1


            sum_lens_p = [ sum_lens_p[i] / counts_p[i] for i in range(0, len(bins))]
            sum_lens_h = [ sum_lens_h[i] / counts_h[i] for i in range(0, len(bins))]


            plt.plot(bins, sum_lens_p)
            plt.xlabel('frquency bin')
            plt.ylabel('mean diff len')
            plt.title('purine. Range: {}-{} bp'.format(range_s, range_e))
            plt.savefig('Great_tit_purine_evolution_range_{}-{}.png'.format(range_s, range_e), dpi=600, bbox_inches='tight')
            plt.close()

            plt.plot(bins, sum_lens_h)
            plt.xlabel('frquency bin')
            plt.ylabel('mean diff len')
            plt.title('hydro. Range: {}-{} bp'.format(range_s, range_e))
            plt.savefig('Great_tit_hydro_evolution_range_{}-{}.png'.format(range_s, range_e), dpi=600, bbox_inches='tight')
            plt.close()

            cur_df = pd.DataFrame({'{}-{}_purine'.format(range_s, range_e) : sum_lens_p, '{}-{}_purine_counts'.format(range_s, range_e) : counts_p,
                                   '{}-{}_purine_counts_g'.format(range_s, range_e): counts_p_g,
                                   '{}-{}_purine_counts_e'.format(range_s, range_e): counts_p_e,
                                   '{}-{}_purine_counts_l'.format(range_s, range_e): counts_p_l,
                                   '{}-{}_hydro'.format(range_s, range_e) : sum_lens_h, '{}-{}_hydro_couns'.format(range_s, range_e) : counts_h,
                                   '{}-{}_hydro_counts_g'.format(range_s, range_e): counts_h_g,
                                   '{}-{}_hydro_counts_e'.format(range_s, range_e): counts_h_e,
                                   '{}-{}_hydro_counts_l'.format(range_s, range_e): counts_h_l
                                   } )
            df_out = pd.concat([df_out, cur_df], axis=1)

    df_out.to_csv('len_effects_Great_tit_binned.tsv', sep='\t', index=False)


def case3():
    #fields = ['variant_id', 'effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard
    #fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard

    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    # fields = ['EAF_A1',  'm1', 'm2', 'p1', 'p2'] #Karlsson
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira

    fields = ['freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit
    #fields = ['CHROM', 'POS', 'freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit


    #df = pd.read_csv('Harlan_hydro.tsv', usecols=fields, sep='\t', header=0)
    df = pd.read_csv('Great_tit_hydro.tsv', usecols=fields, sep='\t', header=0)


    df['m_min'] = [min(x1, x2) for x1, x2 in zip(df.loc[:, 'm1'], df.loc[:, 'm2']) ]
    df['p_min'] = [min(x1, x2) for x1, x2 in zip(df.loc[:, 'p1'], df.loc[:, 'p2']) ]

    perc_purine = int( np.percentile(df.loc[:, 'm_min'], 100) )
    perc_hydro = int( np.percentile(df.loc[:, 'p_min'], 100) )

    counts, bins, bars = plt.hist(df.loc[:, 'm_min'], bins = [x for x in (range(1,  perc_purine) )])
    bins = [ 0.5 * ( bins[i] + bins[i + 1] )  for i in range(len(bins) - 1) ]


    df_out = pd.DataFrame({'counts' : counts, 'bins_center' : bins } )
    df_out.to_csv('hist_Great_tit_purines.tsv', sep='\t', index=False)


    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel('chain length bin')
    plt.ylabel('SNP counts')
    plt.title('Great_tit purine chain length histogram')
    plt.savefig('Great_tit_purine_hist.png', dpi=600, bbox_inches='tight')
    plt.close()

    counts, bins, bars = plt.hist(df.loc[:, 'p_min'], bins = [x for x in (range(1,  perc_hydro) )])
    bins = [ 0.5 * ( bins[i] + bins[i + 1] )  for i in range(len(bins) - 1) ]

    plt.xscale('log')
    plt.yscale('log')

    plt.xlabel('chain length bin')
    plt.ylabel('SNP counts')
    plt.title('Great_tit hydro chain length histogram')
    plt.savefig('Great_tit_hydro_hist.png', dpi=600, bbox_inches='tight')
    plt.close()

    df_out = pd.DataFrame({'counts' : counts, 'bins_center' : bins } )
    df_out.to_csv('hist_Great_tit_hydro.tsv', sep='\t', index=False)


def caseW():
    fields = ['CHROM', 'POS', 'freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit



    df = pd.read_csv('Turcot_hydro_gen.tsv', sep='\t', header=0)

    print(df.iloc[:15, 2:7])

    df = df[df['ExAC_MAF'] != '-']

    print(df.iloc[:15, 2:7])

    def func(row):
        splited0 = str.split(row['ExAC_MAF'], ',')
        found = False

        for k in range(len(splited0)):
            nucl = splited0[k][0]
            after_nucl = splited0[k][1]
            if (nucl == row['REF'] or nucl == row['ALT']) and after_nucl == ':':
                found = True
                found_K  = k
                break

        if not found:
            return np.NaN

        splited = str.split(splited0[found_K], ':')
        allele = splited[0]
        freq = float(splited[1])

        if allele == row['REF']:
            return 1 - freq
        elif allele == row['ALT']:
            return freq
        else:
            print(splited0)
            assert(False)

    df['ExAC_MAF'] = df.apply(func, axis = 1)
    df.dropna(inplace=True, axis=0)

    df.to_csv('Turcot_hydro_gen_reformat.tsv', sep='\t', index=False)

if __name__ == '__main__':
    #caseW()
    #case1()
    #case2()
    #case3()
    #caseAT()
    case_purine_shoulder()
    #case_AT_get()

