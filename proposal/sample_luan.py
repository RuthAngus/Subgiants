import numpy as np

ID = np.genfromtxt('sample_luan.out', skip_header=1, usecols=0, dtype=str).T
mass = np.genfromtxt('sample_luan.out', skip_header=1, usecols=9).T

id_list = ["hd93396", "hd96063", "hd96683", "hd98219", "hd118082", \
        "hd144363", "hd202867", "hd210521", "hd213278", "hd216834"]

masses = []
for i in id_list:
    print i, ID[ID==i][0], mass[ID==i][0]
    masses.append(mass[ID==i][0])

print min(masses), max(masses)
