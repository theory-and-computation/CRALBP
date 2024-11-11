import MDAnalysis as mda
import matplotlib.pyplot as plt
import scienceplots
import matplotlib
import math
import seaborn as sns
import numpy as np
from matplotlib import rc
import pandas as pd


def initializeUniverse(PSF, *DCD) -> mda.Universe:
    """Initializes the universe object."""
    uni = mda.Universe(PSF, DCD)

    return uni


def calcDistance(com1, com2):
    x_diff = abs(com2[0] - com1[0])
    y_diff = abs(com2[1] - com1[1])
    z_diff = abs(com2[2] - com1[2])

    return math.sqrt((x_diff ** 2) + (y_diff ** 2) + (z_diff ** 2))


def getCenterOfMassCoordinates(universe, select: str = "") -> (dict, int):
    """Calculates center of mass coordinates for a protein in a trajectory file."""
    total_frames = len(universe.trajectory)
    selection = universe.select_atoms(f'{select}')

    protein_COM_coordinates = {}

    for frame_index in range(total_frames):

        universe.trajectory[frame_index]
        protein_COM = selection.center_of_mass()

        protein_COM_coordinates[frame_index] = protein_COM

    return protein_COM_coordinates


def getPairWiseDistance(selection1, selection2):

    pairwise_distances = []

    if len(selection1) == len(selection2):
        timesteps = len(selection1)
        for timestep in range(timesteps):
            pairwise_distance = calcDistance(selection1[timestep], selection2[timestep])
            pairwise_distances.append(pairwise_distance)
    else:
        raise Exception("Lists to be compared are not of equal length.")

    return pairwise_distances


def getAllPairwiseDistances(universe=None, pairz=None):
    """
    Formatting data into form:

    distance_dict = {
    "Group 1": [1, 2, 2, 3, 3, 4, 5],
    "Group 2": [2, 3, 3, 4, 4, 5, 6],
    "Group 3": [1, 3, 3, 4, 4, 5, 5],
    "Group 4": [2, 2, 3, 4, 4, 4, 6]
    }
    """
    if not universe:
        return

    all_distances = {}
    for pair in pairz:
        sel1 = getCenterOfMassCoordinates(universe, f"protein and resid {pair[0]}")
        sel2 = getCenterOfMassCoordinates(universe, f"protein and resid {pair[1]}")
        groupLabel = f"{pair[0]}-{pair[1]}"
        p_distances = getPairWiseDistance(sel1, sel2)
        all_distances[groupLabel] = p_distances

    return all_distances


def plotPairWiseDistances(stat=None, tra=None, nm=None, index1=None):
    rc('font', **{'family':'serif', 'serif':['Arial']})
    rc('text', usetex=True)
    fig, ax = plt.subplots()
    ax.set(xlabel="Timesteps", ylabel=f"Distance in Angstroms")
    ax.set(title = 'Pairwise Residue Fluctuations')

    #ex. index1 = "163-215"

    sns.lineplot(x=range(len(nm[index1])), y=nm[index1], color="#b9f2f0", label="Neutral")
    sns.lineplot(x=range(len(tra[index1])), y=tra[index1], color="#d0bbff", label="Transitioned")
    sns.lineplot(x=range(len(stat[index1])), y=stat[index1], color="#ff9f9b", label="Stationary")


    ax.set_xlim([0, 10000])
    figure = ax.get_figure()
    figure.savefig('pairwises_ot.pdf', dpi = 500)
    plt.show()


if __name__ == '__main__':
    matplotlib.rcParams['font.family'] = 'Arial'
    plt.style.use(['science', 'ieee', 'no-latex'])

    psf1 = '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/stripped.psf'
    dcd1 = ['/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00004-stripped.dcd', '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00005-stripped.dcd', '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00006-stripped.dcd', '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00007-stripped.dcd', '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00008-stripped.dcd', '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00009-stripped.dcd']

    psf2 = '/Users/danielsantos/Downloads/CRALBP_analysis/RSMC_stripped_and_aligned/stripped.psf'
    dcd2 = '/Users/danielsantos/Downloads/CRALBP_analysis/RSMC_stripped_and_aligned/rotatedSmallChargedMembrane_2-ions-wats-00004-stripped.dcd'

    psf3 = '/Users/danielsantos/Downloads/CRALBP_analysis/all_stripped_and_aligned/stripped.psf'
    dcd3 = '/Users/danielsantos/Downloads/CRALBP_analysis/all_stripped_and_aligned/rotatedTNeutralMembrane_act-ions-wats-00004-stripped.dcd'

    u1 = initializeUniverse(psf1, dcd1)
    u2 = initializeUniverse(psf2, dcd2)
    u3 = initializeUniverse(psf3, dcd3)

    pairs = [(171, 215), (216, 250), (220, 240), (161, 173), (202, 220), (204, 223), (166, 227), (200, 215), (251, 263), (173, 250), (255, 266), (163, 215)]
    stat_pdist = getAllPairwiseDistances(u1, pairs)
    tra_pdist = getAllPairwiseDistances(u2, pairs)
    nm_pdist = getAllPairwiseDistances(u3, pairs)
    plotPairWiseDistances(stat_pdist, tra_pdist, nm_pdist, index1="163-215")
