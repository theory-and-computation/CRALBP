import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

def readFile(infile_name) -> list:
    with open(infile_name, 'r') as fileobj:
        lines = fileobj.readlines()
        para1 = []
        for line in lines:
            nums = line.split()
            para1.append(int(nums[0]))
        
        return para1

def plot(parameter_set) -> None:
    sns.set()
    fig, ax = plt.subplots(figsize=(15, 8))
    sns.barplot(x=range(len(parameter_set)), y=parameter_set)
    ax.set_xticklabels([])
    plt.show()

if __name__ == '__main__':

    filename = "./R234W_bia.dat"
    parameter = readFile(filename)

    plot(parameter)