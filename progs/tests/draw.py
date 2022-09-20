import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def makePlot():
    file = pd.read_csv("fixed_n_data.csv")
    naive = np.log(file['naive'])
    kmp = np.log(file['kmp'])

    x_axis = np.log([16,32,64,128,256,512,1024,2048,4096,8192,16384,32768])
    plt.scatter(x_axis ,naive, c="blue", label="naive")
    plt.scatter(x_axis ,kmp, c="green", label="kmp")
    plt.xlabel("length of x")
    plt.ylabel("time")
    
    plt.show()



makePlot()