#include <cstdio>

extern "C"
void calc_corr_dna(double* seq, int n_seq, double* corr_func, int n_corr_func)
{
    printf("calc_corr_dna in cpp module\n");
    printf("%d %d \n", n_seq, n_corr_func);

    int test = 1.000000 * 1.000000;

    printf("%d\n", test);

    int was_gap = 0;

    for (int i = 0; i < n_corr_func; ++i)
    {
        int total_sum = 0;
        int total_N = 0;

        for (int k = 0; k < n_corr_func; ++k)
        {
            if (k + i >= n_seq) {
                break;
            }

            if ( (seq[k] == -1) || (seq[k + i] == -1) ){   // -1 is gap indicator
                was_gap = 1;
                continue;
            }

            total_N++;

//            total_sum += seq[k] * seq[k + i];

            if (seq[k] == seq[k + i]){
               total_sum++;
            }
        }

        if (total_N > 0){
            corr_func[i] = (double (total_sum) ) / total_N;
        }
    }

    printf("%d \n", was_gap);

}
