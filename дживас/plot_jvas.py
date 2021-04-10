import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    fields = ['other_allele', 'effect_allele', 'effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    # fields = ['EAF_A1',  'm1', 'm2', 'p1', 'p2'] #Karlsson
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira
    #fields = ['A1', 'A2', 'freq2',  'm1', 'm2', 'p1', 'p2'] #freq2, GreatTit


    df = pd.read_csv('Howard_test_hydro_ext2.tsv', usecols=fields, sep='\t', header=0)

    filters = [['A', 'G'], ['A','C'], ['A', 'T'], ['C','G'], ['C','T'], ['G', 'T']]

    #filters = [['A', 'G', 'C', 'T']]
    print(df.head)
    for min_diff in [0, 2, 4, 6]:
        for filter in filters:
            points = np.arange(0.5, 1.01, 0.05)
            bins = [ 0.5 * ( points[i] + points[i + 1] ) for i in range(0, len(points) - 1 )]

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
            plt.savefig('Howard_purine_evolution_diff_{}_filt_{}.png'.format(min_diff, filter[0] + filter[1]), dpi=600, bbox_inches='tight')
            plt.close()

            plt.plot(bins, sum_lens_h)
            plt.xlabel('frquency bin')
            plt.ylabel('mean diff len')
            plt.title('hydro. Min Diff: {} bp'.format(min_diff))
            plt.savefig('Howard_hydro_evolution_diff_{}_filt_{}.png'.format(min_diff, filter[0] + filter[1]), dpi=600, bbox_inches='tight')
            plt.close()

        #df_out = pd.DataFrame( {'bins' : bins, 'counts' :  counts, 'mean_diff_p' : sum_lens_p, 'mean_diff_h' : sum_lens_h} )
        #df_out.to_csv('effect.tsv', sep='\t', index=False)