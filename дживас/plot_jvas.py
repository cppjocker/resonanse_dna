import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Howard
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #deLang
    # fields = ['EAF_A1',  'm1', 'm2', 'p1', 'p2'] #Karlsson
    # fields = ['effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Wojcik
    # fields = ['hm_effect_allele_frequency',  'm1', 'm2', 'p1', 'p2'] #Fereira
    # fields = ['freq2',  'm1', 'm2', 'p1', 'p2'] #freq2


    df = pd.read_csv('Howard_test_chunks.tsv', usecols=fields, sep='\t', header=0)
#    print(df.head)
    points = np.arange(0.5, 1.01, 0.05)
    bins = [ 0.5 * ( points[i] + points[i + 1] ) for i in range(0, len(points) - 1 )]

    counts = np.zeros( len(bins) )
    sum_lens_p = np.zeros( len(bins) )
    sum_lens_h = np.zeros( len(bins) )

    for k in range(0, df.shape[0]):
        freq2 = df.iloc[k, 0]
        freq1 = 1 - freq2

        m1 = df.iloc[k, 1]
        m2 = df.iloc[k, 2]

        p1 = df.iloc[k, 3]
        p2 = df.iloc[k, 4]

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
                counts[i] += 1
                sum_lens_p[i] += diff_m
                sum_lens_h[i] += diff_h

    print(bins)
    print(counts)
    print(sum_lens_p)
    print(sum_lens_h)


    sum_lens_p = [ sum_lens_p[i] / counts[i] for i in range(0, len(bins))]
    sum_lens_h = [ sum_lens_h[i] / counts[i] for i in range(0, len(bins))]


    plt.plot(bins, sum_lens_p)
    plt.xlabel('frquency bin')
    plt.ylabel('mean diff len')
    plt.title('purine')
    plt.savefig('Howard_test_purine_evolution.png', dpi=600, bbox_inches='tight')
    plt.close()

    plt.plot(bins, sum_lens_h)
    plt.xlabel('frquency bin')
    plt.ylabel('mean diff len')
    plt.title('hydro')
    plt.savefig('Howard_test_hydro_evolution.png', dpi=600, bbox_inches='tight')
    plt.close()

    df_out = pd.DataFrame( {'bins' : bins, 'counts' :  counts, 'mean_diff_p' : sum_lens_p, 'mean_diff_h' : sum_lens_h} )
    df_out.to_csv('effect.tsv', sep='\t', index=False)