code_list = ['hd_07', 'hd_71', 'hd_72', 'hd_73', 'hd_74', 'hd_75', 'hd_76', 'hd_77', 'hd_78', 'hd_79', 'hd_80', 'hd_81', 'hd_82', 'hd_83', 'hd_84', 'hd_85', 'hd_86']

code_list = ['hd_07', 'hd_55', 'hd_56', 'hd_57', 'hd_58', 'hd_59', 'hd_60', 'hd_61', 'hd_62', 'hd_63', 'hd_64', 'hd_65', 'hd_66', 'hd_67', 'hd_68', 'hd_69', 'hd_70']


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # fields = ['effect_allele_frequency'] #Howard
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    fields = ['EAF_A1'] #Karlsson
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira
    # fields = ['freq2',  'm1', 'm2', 'p1', 'p2'] #freq2

    fields += code_list
    fields += ['purine']

    df = pd.read_csv('Karlsson_new_many_hydro.tsv', usecols=fields, sep='\t', header=0)


#    print(df.head)
    points = np.arange(0.5, 1.01, 0.05)
    bins = [ 0.5 * ( points[i] + points[i + 1] ) for i in range(0, len(points) - 1 )]
    df_out = pd.DataFrame( {'bins' : bins} )


    for range_len in [ [1, 1], [2, 2], [3, 3], [4, 5], [6, 10], [11, 15], [16, 10000] ]:
    #for range_len in [[1, 10000]]:
        range_s = range_len[0]
        range_e = range_len[1]

        counts = np.zeros( len(bins) )
        sum_lens_p = np.zeros( len(bins) )
        sum_lens_h = np.zeros( (len(bins), len(code_list)) )

        counts_h = np.zeros( (len(bins), len(code_list)))
        counts_h_g = np.zeros( (len(bins), len(code_list)))
        counts_h_e = np.zeros( (len(bins), len(code_list)))
        counts_h_l = np.zeros( (len(bins), len(code_list)))


        for k in range(0, df.shape[0]):
            freq2 = df.iloc[k, 0]
            freq1 = 1 - freq2

            purine = df.iloc[k, -1]
            m1 = purine.split(';')[0]
            m2 = purine.split(';')[1]

            m1 = int(m1)
            m2 = int(m2)

            freq = max(freq1, freq2)
            diff_h = np.zeros(len(code_list))
            min_p = np.zeros(len(code_list))

            if freq1 > freq2:
                freq = freq1
                diff_m = m1 - m2
                for j in range(0, len(code_list)):
                    hydro = df.iloc[k, 1 + j]
                    p1 = hydro.split(';')[0]
                    p2 = hydro.split(';')[1]

                    p1 = int(p1)
                    p2 = int(p2)

                    diff_h[j] = p1 - p2
                    min_p[j] = min(p1, p2)

            else:
                freq = freq2
                diff_m = m2 - m1
                for j in range(0, len(code_list)):
                    hydro = df.iloc[k, 1 + j]
                    p1 = hydro.split(';')[0]
                    p2 = hydro.split(';')[1]

                    p1 = int(p1)
                    p2 = int(p2)

                    diff_h[j] = p2 - p1
                    min_p[j] = min(p2, p1)

            for i in range(0, len(bins) ):
                if( points[i] <= freq  and freq < points[i + 1] ):
                    counts[i] += 1
                    sum_lens_p[i] += diff_m
                    for j in range(0, len(code_list)):
                        if range_s <= min_p[j] and min_p[j] <= range_e:
                            sum_lens_h[i, j] += diff_h[j]
                            counts_h[i, j] += 1

                            if diff_h[j] > 0:
                                counts_h_g[i, j] += 1
                            elif diff_h[j] == 0:
                                counts_h_e[i, j] += 1
                            else:
                                counts_h_l[i, j] += 1


        #print(bins)
        #print(counts)
        #print(sum_lens_p)
        #print(sum_lens_h)


        sum_lens_p = [ sum_lens_p[i] / counts[i] for i in range(0, len(bins))]


        plt.plot(bins, sum_lens_p)
        plt.xlabel('frequency bin')
        plt.ylabel('mean diff len')
        plt.title('Purine. Range {}-{}'.format(range_s, range_e))
        plt.savefig('Karlsson_new_many_purine_evolution_{}-{}.png'.format(range_s, range_e), dpi=600, bbox_inches='tight')
        plt.close()

        markers = ['*', '^', 'o', '+', 'v']
        for j in range(0, len(code_list)):
            sum_lens_h[:, j] = [ sum_lens_h[i, j] / counts[i] for i in range(0, len(bins))]

            marker = markers[j % 5]
            plt.plot(bins, sum_lens_h[:, j], marker=marker)

            cur_df = pd.DataFrame(data=sum_lens_h[:, j], columns=[code_list[j]], index=df_out.index)
            cur_df = pd.DataFrame(
                {                 '{}-{}_hydro_{}'.format(range_s, range_e, code_list[j]): sum_lens_h[:, j],
                 '{}-{}_hydro_{}_counts'.format(range_s, range_e, code_list[j]): counts_h[:, j],
                 '{}-{}_hydro_{}_counts_g'.format(range_s, range_e, code_list[j]): counts_h_g[:, j],
                 '{}-{}_hydro_{}_counts_e'.format(range_s, range_e, code_list[j]): counts_h_e[:, j],
                 '{}-{}_hydro_{}_counts_l'.format(range_s, range_e, code_list[j]): counts_h_l[:, j]
                 })
            df_out = pd.concat( [df_out, cur_df  ],  axis=1)

        plt.legend(code_list, bbox_to_anchor=(1.05, 1),  loc='upper left')
        plt.xlabel('frequency bin')
        plt.ylabel('mean diff len')
        plt.title('Hydro. Range {}-{}'.format(range_s, range_e))
        plt.savefig('Karlsson_hydro_new_many_evolution_{}-{}.png'.format(range_s, range_e), dpi=600, bbox_inches='tight')
        plt.close()



    df_out.to_csv('Karlsson_new_many_hydro_evolution_binned.tsv', sep='\t', index=False)