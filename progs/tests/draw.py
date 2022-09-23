import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def makePlot():
    
    file = pd.read_csv("same_alg_twice.csv")
    naive = np.log(file['kmp1'])
    kmp = np.log(file['kmp2'])
    labels = ['naive', 'kmp']

    x_axis = np.log((file['x_size']))
    plt.scatter(x_axis ,naive, c="blue", label="naive, p=100")
    plt.scatter(x_axis ,kmp, c="green", label="naive, p=200")
    plt.legend(loc="upper left")
    plt.xlabel(" log (length of x)")
    plt.ylabel("log (time)")
    
    #plt.ticklabel_format(useOffset=False, style='plain')
    
    plt.show()



makePlot()