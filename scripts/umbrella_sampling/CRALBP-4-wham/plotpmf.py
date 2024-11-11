import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rc
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import math
import scienceplots

def readFile(infile_name) -> (list, list):
    with open(infile_name, 'r') as fileobj:
        lines = fileobj.readlines()
        para1 = []
        para2 = []
        for line in lines:
            nums = line.split()
            para1.append(float(nums[0]))
            para2.append(float(nums[1]))
        
        return para1, para2


def plot(pmf1) -> None:
    plt.style.use(['science', 'ieee', 'muted'])
    rc('font',**{'family':'serif','serif':['Arial']})
    rc('text', usetex=True)
    fig, ax = plt.subplots(figsize=(15, 8))

    plt.plot(pmf1[0], pmf1[1], '-', label='S')


    ax.set_xlim(-10, 10)

    figure = ax.get_figure()
    figure.savefig('R234W_pmf.dat.pdf')

    plt.legend()
    plt.show()

if __name__ == '__main__':

    filename1 = './R234W_pmf.dat'
    parameters1 = readFile(filename1)
    plot(parameters1)
    