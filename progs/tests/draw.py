import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def makePlot():
    
    file = pd.read_csv("fixed_n_data.csv")
    naive = (file['naive'])
    kmp = file['kmp']
    labels = ['naive', 'kmp']

    x_axis = ((file['x_size']))
    #plt.scatter(x_axis ,naive, c="blue", label="naive")
    plt.scatter(x_axis ,kmp, c="green", label="kmp")
    plt.legend(loc="upper left")
    plt.xlabel("length of x")
    plt.ylabel("time")
    
    #plt.ticklabel_format(useOffset=False, style='plain')
    
    plt.show()



makePlot()