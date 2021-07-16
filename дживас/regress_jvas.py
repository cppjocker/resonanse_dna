import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import seaborn as sns

def case1():
    use_counts = True

    inc = 1
    if use_counts:
        inc = 2

    df = pd.read_csv('len_effects_great_tit.tsv', sep='\t', header=0)
    start_col = 1

    bins = df.iloc[:, 0]
    X = bins[:, None]

    angles = []
    r2s = []
    total_counts = []
    h_columns = []
    for i in range(start_col, df.shape[1], inc):
        cur_column = df.columns[i]
        cur_h = df.iloc[:, i].values

        if( any( np.isnan( cur_h  ) ) ):
            continue


        reg = LinearRegression()
        reg .fit(X, cur_h)
        y_pred = reg.predict(X)
        r2 = r2_score(cur_h, y_pred)

        slope = reg.coef_[0]

        angle = (math.atan(slope)) / math.pi * 180
        angles.append(angle)
        r2s.append(r2 * 100)

        h_columns.append(cur_column)

        if(use_counts):
            cur_counts = sum(df.iloc[:, i + 1].values)
            total_counts.append(cur_counts)

    angles1 = pd.DataFrame({'option' : h_columns, 'value' : angles, 'type' : ['angle'] * len(h_columns)})
    angles2 = pd.DataFrame({'option' : h_columns, 'value' : r2s, 'type' : ['r2'] * len(h_columns)})
    angles3 = pd.DataFrame({'option' : h_columns, 'value' : total_counts, 'type' : ['counts'] * len(h_columns)})

    angles = pd.concat([angles1, angles2, angles3], axis=0)

    ax = sns.barplot(x="option", y="value", data=angles, hue="type")

    for item in ax.get_xticklabels():
        item.set_rotation(45)

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=3)

    plt.savefig("len_Great_Tit_angles.png", bbox_inches = 'tight', dpi = 600)
    plt.close()

    angles.to_csv("len_Great_Tit_angles.tsv", sep = '\t', index=False)


def case2():

    inc = 10

    df = pd.read_csv('len_effects_Howard_hd87_binned.tsv', sep='\t', header=0)
    fontsize_x = 2

    bins = df.iloc[:, 0]


    chain_len = 0

    cdiffs = pd.DataFrame({'bin_chain_len' : [], 'count_diff_len' : [], 'type' : []})

    purine_graph = False

    if purine_graph:
        desc = 'Purine\n\n'
        start_col = 1
        file = 'len_effects_Howard_purine_binned_diff'
    else:
        start_col = 6
        desc = 'Hydro\n\n'
        file = 'len_effects_Howard_hydro_binned_diff'

    empty_str = ' '
    empty_df = pd.DataFrame({'bin_chain_len': [empty_str], 'count_diff_len': [0.1], 'type': '{:.2f}'.format(bins[0])})

    cdiffs = pd.concat([cdiffs, empty_df])

    for i in range(start_col, df.shape[1], inc):
        col = str.split(df.columns[i], '_')[0]
        cur_desc = '{0} id --- {1} Range bp\n'.format(chain_len, col)
        desc += cur_desc

        count_g = df.iloc[:, i + 2]
        count_l = df.iloc[:, i + 4]
        count_diff = count_g - count_l


        if( any( np.isnan( count_diff  ) )   ):
            assert(False)

        chain_id = [ '{0}_{1}'.format(chain_len, j) for j in range(len(bins))]
        bins_id = [ '{:.2f}'.format(bins[j]) for j in range(len(bins))]

        cur_diff = pd.DataFrame({'bin_chain_len' : chain_id, 'count_diff_len' : count_diff, 'type' : bins_id})

        empty_str += ' '
        empty_df = pd.DataFrame({'bin_chain_len': [empty_str], 'count_diff_len': [0.1], 'type': '{:.2f}'.format(bins[0])})

        cdiffs = pd.concat([cdiffs,  cur_diff, empty_df] )

        chain_len+=1

    empty_str += ' '
    empty_df = pd.DataFrame({'bin_chain_len': [empty_str], 'count_diff_len': [0.1], 'type': '{:.2f}'.format(bins[0])})


    cdiffs = pd.concat([cdiffs, empty_df])

    cdiffs.reset_index(inplace=True)
    ax = sns.barplot(x="bin_chain_len", y="count_diff_len", data=cdiffs, hue = 'type')

    cdiffs = cdiffs.loc[cdiffs['count_diff_len'] != 0.1, :]

    ax.legend(bbox_to_anchor=(1, 1), loc=2)
    ax.set_title(desc)

    fig = plt.gcf()
    fig.canvas.draw()

    for item in ax.get_xticklabels():
       item.set_rotation(45)

    ax.set_xticklabels(ax.get_xticklabels(), fontsize=fontsize_x)
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=5)

    plt.savefig( file + ".png", bbox_inches = 'tight', dpi = 600)
    plt.close()

    cdiffs.to_csv( file + ".tsv", sep = '\t', index=False)


def case3():
    files = ["len_effects_Howard_merged.tsv", "len_effects_Karlsson_merged.tsv",
             "len_effects_Wojcik_merged.tsv", "len_effects_Ferreira_merged.tsv", "len_effects_LeeJJ_merged.tsv"]

    names = ["Howard", "Karlsson", "Wojcik", "Ferreira", "LeeJJ"]

    for i in range(len(files) ):
        file = files[i]
        df = pd.read_csv(file, sep='\t', header=0)
        x = df.iloc[:, 0]
        y = df.iloc[:, 1]
        plt.plot(x, y, label = names[i])

        if i == 0:
            df_out = pd.DataFrame(columns=['bins'], data=x.values)

        df_out_cur = pd.DataFrame(columns=[names[i]], data=y.values)

        df_out = pd.concat( [df_out, df_out_cur], axis = 1 )

    plt.xlabel('frequency')
    # Set the y axis label of the current axis.
    plt.ylabel('average length diff')
    # Set a title of the current axes.
    plt.title('Purine evolution')
    # show a legend on the plot
    plt.legend()

    plt.savefig("human_purine_many.png", bbox_inches='tight', dpi=600)
    plt.close()

    df_out.to_csv("human_purine_many.tsv", sep='\t', index=False)

    for i in range(len(files) ):
        file = files[i]
        df = pd.read_csv(file, sep='\t', header=0)
        x = df.iloc[:, 0]
        y = df.iloc[:, 6]
        plt.plot(x, y, label = names[i])

        if i == 0:
            df_out = pd.DataFrame(columns=['bins'], data=x.values)

        df_out_cur = pd.DataFrame(columns=[names[i]], data=y.values)

        df_out = pd.concat( [df_out, df_out_cur], axis = 1 )


    plt.xlabel('frequency')
    # Set the y axis label of the current axis.
    plt.ylabel('average length diff')
    # Set a title of the current axes.
    plt.title('Hydro evolution')
    # show a legend on the plot
    plt.legend()

    plt.savefig("human_hydro_many.png", bbox_inches='tight', dpi=600)
    plt.close()

    df_out.to_csv("human_hydro_many.tsv", sep='\t', index=False)


def case4():
    files = ["len_effects_Howard_binned.tsv", "len_effects_Karlsson_binned.tsv",
             "len_effects_Wojcik_binned.tsv", "len_effects_Ferreira_binned.tsv", "len_effects_LeeJJ_binned.tsv",
             "len_effects_Mahajan_binned.tsv", "len_effects_Turcot_binned.tsv",
             "len_effects_salmon_binned.tsv", "len_effects_Great_tit_binned.tsv"]

    names = ["Howard", "Karlsson", "Wojcik", "Ferreira", "LeeJJ", "Mahajan", "Turcot", "salmon", "Great_tit"]


    purine_graph = True
    inc = 10
    empty_str = ' '


    if purine_graph:
        desc = 'Purine\n\n'
        start_col = 1
    else:
        start_col = 6
        desc = 'Hydro\n\n'


    all_diff_species = pd.DataFrame({'bin_chain_len': [], 'count_diff_len': [], 'type' : []})

    for i in range(len(files) ):
        file = files[i]
        df = pd.read_csv(file, sep='\t', header=0)
        bins = df.iloc[:, 0]
        chain_len = 0

        # if i == 0:
        #     df_out = pd.DataFrame(columns=['bins'], data=x.values)
        #
        # df_out_cur = pd.DataFrame(columns=[names[i]], data=y.values)
        #
        # df_out = pd.concat( [df_out, df_out_cur], axis = 1 )

        all_diff = pd.DataFrame({'bin_chain_len': [], 'count_diff_len': []})

        for j in range(start_col, df.shape[1], inc):
            col = str.split(df.columns[j], '_')[0]
            cur_desc = '{0} id --- {1} Range bp\n'.format(chain_len, col)

            if i == 0:
                desc += cur_desc

            count_g = df.iloc[:, j + 2]
            count_l = df.iloc[:, j + 4]
            count_diff = count_g - count_l

            if (any(np.isnan(count_diff))):
                assert (False)

            chain_id = ['{0}_{1}'.format(chain_len, j) for j in range(len(bins))]

            cur_diff = pd.DataFrame({'bin_chain_len': chain_id, 'count_diff_len': count_diff})

            all_diff = pd.concat([all_diff, cur_diff], axis=0)


            chain_len += 1

        plt.plot(all_diff.loc[:, 'bin_chain_len'], all_diff.loc[:, 'count_diff_len'], label=names[i])
        diff_species = pd.DataFrame({'bin_chain_len': all_diff.loc[:, 'bin_chain_len'],
                                     'count_diff_len': all_diff.loc[:, 'count_diff_len'],
                                                       'type' : [names[i]] * all_diff.shape[0]})

        all_diff_species = pd.concat([all_diff_species, diff_species], axis=0)

    plt.xlabel('frequency')
    # Set the y axis label of the current axis.
    plt.ylabel('counts splined')
    # Set a title of the current axes.
    plt.title('Purine evolution')
    # show a legend on the plot
#    plt.legend()

    ax = plt.gca()
    ax.legend(bbox_to_anchor=(1, 1), loc=2)


    ax.set_title(desc)
    ax.set_yscale('symlog')

    ax.tick_params(axis="x", labelsize=2)

    #ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
    #ax.set_yticklabels(ax.get_yticklabels(), fontsize=25)


    plt.savefig("human_purine_many_counts.png", bbox_inches='tight', dpi=600)
    plt.close()

    all_diff_species.to_csv("human_purine_many_counts.tsv", sep='\t', index=False)



#Many Hydro
def case5():
    file = "Great_Tit_new_many_hydro_evolution_binned.tsv"


    start_col = 1
    desc = 'Hydro\n\n'

    all_diff_species = pd.DataFrame({'bin_chain_len': [], 'count_diff_len': [], 'type' : []})

    df = pd.read_csv(file, sep='\t', header=0)

    code_list = ['hd_07', 'hd_55', 'hd_56', 'hd_57', 'hd_58', 'hd_59', 'hd_60', 'hd_61', 'hd_62', 'hd_63', 'hd_64',
                 'hd_65', 'hd_66', 'hd_67', 'hd_68', 'hd_69', 'hd_70']

    inc = 5 * len(code_list)


    df = pd.read_csv(file, sep='\t', header=0)
    bins = df.iloc[:, 0]

    for i in range(len(code_list) ):
        chain_len = 0

        all_diff = pd.DataFrame({'bin_chain_len': [], 'count_diff_len': []})

        for j in range(start_col, df.shape[1], inc):
            col = str.split(df.columns[j], '_')[0]
            cur_desc = '{0} id --- {1} Range bp\n'.format(chain_len, col)

            if i == 0:
                desc += cur_desc

            count_g = df.iloc[:, j + 2]
            count_l = df.iloc[:, j + 4]
            count_diff = count_g - count_l

            print(j)

            if (any(np.isnan(count_diff))):
                assert (False)

            chain_id = ['{0}_{1}'.format(chain_len, j) for j in range(len(bins))]

            cur_diff = pd.DataFrame({'bin_chain_len': chain_id, 'count_diff_len': count_diff})

            all_diff = pd.concat([all_diff, cur_diff], axis=0)


            chain_len += 1

        plt.plot(all_diff.loc[:, 'bin_chain_len'], all_diff.loc[:, 'count_diff_len'], label=code_list[i])
        diff_species = pd.DataFrame({'bin_chain_len': all_diff.loc[:, 'bin_chain_len'],
                                     'count_diff_len': all_diff.loc[:, 'count_diff_len'],
                                                       'type' : [code_list[i]] * all_diff.shape[0]})

        all_diff_species = pd.concat([all_diff_species, diff_species], axis=0)

        start_col += 5

    plt.xlabel('frequency')
    # Set the y axis label of the current axis.
    plt.ylabel('counts splined')
    # Set a title of the current axes.
    plt.title('Hydro evolution. Salmon')
    # show a legend on the plot
#    plt.legend()

    ax = plt.gca()
    ax.legend(bbox_to_anchor=(1, 1), loc=2)


    ax.set_title(desc)
    ax.set_yscale('symlog')

    ax.tick_params(axis="x", labelsize=2)

    #ax.set_xticklabels(ax.get_xticklabels(), fontsize=15)
    #ax.set_yticklabels(ax.get_yticklabels(), fontsize=25)


    plt.savefig("Great_Tit_many_hydro_binned.png", bbox_inches='tight', dpi=600)
    plt.close()

    all_diff_species.to_csv("Great_Tit_many_hydro_binned.tsv", sep='\t', index=False)



if __name__ == '__main__':
    #case1()
    #case2()
    #case3()
    #case4()
    case5()