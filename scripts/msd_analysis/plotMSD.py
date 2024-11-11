import numpy as npy
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')
def openMSDFile(filename):

    with open(filename, 'rb') as file:

        datasets = npy.load(file)

    return datasets


def plotMSD(MSD_DATA):
    """Create plots of the data."""

    chargedMembraneMSD = MSD_DATA[0]
    CM_stdev = MSD_DATA[2]

    neutralMembraneMSD = MSD_DATA[1]
    NM_stdev = MSD_DATA[3]

    tau_to_ns = [x * 20000 * 0.000001 for x in range(len(NM_stdev))]

    popt, pcov = curve_fit(funcToFit, tau_to_ns, neutralMembraneMSD, p0=[10, 0.5], maxfev = 10000)
    popt2, pcov2 = curve_fit(funcToFit, tau_to_ns, chargedMembraneMSD, p0=[10, 2], maxfev = 10000)
    testing_y = [funcToFit(x, *popt) for x in tau_to_ns]
    testing_y2 = [funcToFit(x, *popt2) for x in tau_to_ns]
    print(popt)
    print(popt2)

    NM_upper_stdev = [((NM_stdev[x] / 4) + testing_y[x]) for x in range(len(NM_stdev))]
    NM_lower_stdev = [(testing_y[x] - (NM_stdev[x] / 4)) for x in range(len(NM_stdev))]
    CM_upper_stdev = [((CM_stdev[x] / 4) + testing_y2[x]) for x in range(len(CM_stdev))]
    CM_lower_stdev = [(testing_y2[x] - (CM_stdev[x] / 4)) for x in range(len(CM_stdev))]

    sns.set(font_scale=1.5)
    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(15, 8))

    ax.spines['left'].set_linewidth(2)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['bottom'].set_color('black')

    jump = 850
    ax.fill_between(x=tau_to_ns, y1=CM_upper_stdev, y2=CM_lower_stdev, color='green', alpha=0.2)
    ax.fill_between(x=tau_to_ns, y1 = NM_upper_stdev, y2 = NM_lower_stdev, color = 'green', alpha = 0.2)
    ax.scatter(x=tau_to_ns[::jump], y=chargedMembraneMSD[::jump], s=125, facecolors='none', edgecolors='red', label='POPS species absent')
    ax.scatter(x=tau_to_ns[::jump], y=neutralMembraneMSD[::jump], s=125, facecolors='none', edgecolors='green', label='POPS species present')

    sns.lineplot(x=tau_to_ns, y=testing_y, label='POPS species present - fit')
    sns.lineplot(x=tau_to_ns, y=testing_y2, label='POPS species absent - fit')

    #ax.set(xlabel="ns", ylabel=f"MSD <r²({tau})>")
    #ax.tick_params(axis='y', labelsize = 12, labelleft=False)
    #ax.tick_params(axis = 'x', labelsize = 15, pad=10, labelbottom=False)
    ax.tick_params(bottom = True, top = True, left = True, right = True)
    ax.tick_params(axis='y', labelsize = 12)
    ax.tick_params(axis = 'x', labelsize = 15, pad=10)

    plt.legend().remove()

    figure = ax.get_figure()
    figure.savefig('MSD_test.pdf',bbox_inches='tight', dpi = 500)

    plt.show()


def funcToFit(tau, diff_coef, alpha):
    return diff_coef * (tau ** alpha)


if __name__ == '__main__':

    _data = openMSDFile('/Users/danielsantos/Downloads/CRALBP/save_files/MSD.npy')
    plotMSD(_data)