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
        
        pair = [para1, para2]
        return pair


def calculateErrorCorridor(_mean, observed):
    observed = np.array(observed)
    _mean = np.array(_mean)

    error = (observed - _mean)**2
    half_error = [x/2 for x in error]
    error_corridor_high = observed + half_error
    error_corridor_low = observed - half_error

    return error_corridor_high, error_corridor_low


def calculateVariance(pmf1, pmf2, pmf3, pmf4, pmf5):
    per_point_variance = []
    for i in range(len(pmf1[0])):
        block_points = [pmf1[1][i], pmf2[1][i], pmf3[1][i], pmf4[1][i], pmf5[1][i]]
        variance = np.std(block_points)
        print(variance)
        per_point_variance.append(variance)

    return np.array(per_point_variance)


def setWindow(pmf, start, stop):
    for i in range(len(pmf[0])):
        if pmf[0][i] <= start or pmf[0][i] > stop:
            pmf[0][i] = 999
            pmf[1][i] = 999
    pmf[0] = [x for x in pmf[0] if x != 999]
    pmf[1] = [x for x in pmf[1] if x != 999]


def setWindows(pmfs, start, stop):
    for pmf in pmfs:
        setWindow(pmf, start, stop)

def plotBlocks(pmf1, pmf2, pmf3, pmf4, pmf5, true_pmf):
    plt.style.use(['science', 'ieee', 'muted'])
    rc('font',**{'family':'serif','serif':['Arial']})
    rc('text', usetex=True)
    fig, ax = plt.subplots(figsize=(15, 8))

    setWindows([pmf1, pmf2, pmf3, pmf4, pmf5, true_pmf], start=-10, stop=10)
    total_variance = calculateVariance(pmf1, pmf2, pmf3, pmf4, pmf5)
    half_variance = total_variance / 2

    plt.plot(true_pmf[0], true_pmf[1], '-', label='Pore 1 - POPS bound')
    plt.fill_between(true_pmf[0], true_pmf[1] - half_variance, true_pmf[1] + half_variance, alpha=0.1)

    ax.set_ylim(0, 115)
    ax.set_xlim(-20, 0)
    ax.invert_xaxis()
    plt.legend()
    figure = ax.get_figure()
    figure.savefig('R234_LPUB.pdf')
    plt.show()



def calculateMeanPMF(block1, block2, block3, block4, block5):
    mean_pmf = []
    for i in range(len(block1[1])):
        mean_pmf.append(np.mean([block1[1][i],block2[1][i],block3[1][i],block4[1][i],block5[1][i]]))

    return mean_pmf


if __name__ == '__main__':

    plt.style.use(['science', 'ieee', 'muted'])

    filename1 = './234lpub_block1_pmf.dat'
    filename2 = './234lpub_block2_pmf.dat'
    filename3 = './234lpub_block3_pmf.dat'
    filename4 = './234lpub_block4_pmf.dat'
    filename5 = './234lpub_block5_pmf.dat'
    filename6 = './234lpub_pmf.dat'


    parameters1 = readFile(filename1)
    parameters2 = readFile(filename2)
    parameters3 = readFile(filename3)
    parameters4 = readFile(filename4)
    parameters5 = readFile(filename5)
    parameters6 = readFile(filename6)


    plotBlocks(parameters1, parameters2, parameters3, parameters4, parameters5, parameters6)
