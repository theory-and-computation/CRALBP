import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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

def formatDataFrame(parameter_sets) -> pd.DataFrame:
    df = pd.DataFrame(parameter_sets, columns=['Parameter 1, Parameter 2'])
    return df

def plot(parameter_sets) -> None:
    sns.set()
    fig, ax = plt.subplots(figsize=(15, 8))

    sns.lineplot(x=parameter_sets[0], y=parameter_sets[1])
    plt.show()

if __name__ == '__main__':

    filename = "./R234W_rho.dat"
    parameters = readFile(filename)

    plot(parameters)