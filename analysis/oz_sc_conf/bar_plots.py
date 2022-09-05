import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

DIR = "/mnt/perkinsdata/tongue_STOmics/benchmarking/25_2k_4"
FILE = "accuracy.csv"


MEASURES = [("tf", "Recall"), ("tf", "Precision"), ("ft", "Entropy_wo")]
ALGOS = ["Watershed", "Cellpose", "DeepCell"]

# load data
data = pd.read_csv(os.path.join(DIR, FILE))

# get plotting data from data
means = []
stds = []
for m in MEASURES:
    means.append(list(data[data["TF_FT"] == m[0]][m[1] + "_m"]))
    stds.append(list(data[data["TF_FT"] == m[0]][m[1] + "_s"]))

X = np.arange(3)

fig = plt.figure()
ax = fig.add_axes([0,0,1,1])
ax.bar(X + 0.00, data[0], color = 'b', width = 0.25)
ax.bar(X + 0.25, data[1], color = 'g', width = 0.25)
ax.bar(X + 0.50, data[2], color = 'r', width = 0.25)