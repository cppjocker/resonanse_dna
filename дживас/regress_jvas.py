import pandas as pd
import math
import numpy as np
import matplotlib.pyplot as plt

from sklearn.linear_model import LinearRegression
import seaborn as sns

df = pd.read_csv('effect_many_h.tsv', sep='\t', header=0)

bins = df.iloc[:, 0]
X = bins[:, None]

angles = []
h_columns = []
for i in range(3, df.shape[1]):
    cur_column = df.columns[i]
    cur_h = df.iloc[:, i].values


    reg = LinearRegression().fit(X, cur_h)
    slope = reg.coef_[0]

    angle = (math.atan(slope)) / math.pi * 180
    angles.append(angle)
    h_columns.append(cur_column)


hydro_angles = pd.DataFrame({'code' : h_columns, 'angle' : angles})

ax = sns.barplot(x="code", y="angle", data=hydro_angles, color="blue")

for item in ax.get_xticklabels():
    item.set_rotation(45)

plt.savefig("hydro_angles.png", bbox_inches = 'tight', dpi = 600)
plt.close()