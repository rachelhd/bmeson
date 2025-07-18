import numpy as np
import pandas as pd
from csv import writer


file_dir='/home/rachel/PhD/Bmeson/hl_corrs'
#with open(file_dir+'/hl_ps_16x128_new2.csv') as file:
data = pd.read_csv(file_dir+'/hl_ps_24x40_s0.csv',  header=None)
#    data = data.to_numpy()
#file.close()
"""
data=data[:,1:]
uniquerow = []
repeatedrow = []

for row in data:
    if (row in uniquerow):
        repeatedrow.append(row)
    else: 
        uniquerow.append(row)

print(len(repeatedrow[:,0]))
"""
data_unique = data.drop_duplicates(subset=data.columns[1:])

data_unique.to_csv(file_dir+"/hl_ps_24x40_new2.csv", index=False)
"""
with open(file_dir+"/hl_ps_16x128_new3.csv", 'a') as file:
    write = writer(file)
    for row in uniquerow:
        write.writerow(row)

    file.close
"""

