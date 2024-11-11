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


def getAllPairwiseDistances(universe, pairz):
    """
    Formatting data into form:

    distance_dict = {
    "Group 1": [1, 2, 2, 3, 3, 4, 5],
    "Group 2": [2, 3, 3, 4, 4, 5, 6],
    "Group 3": [1, 3, 3, 4, 4, 5, 5],
    "Group 4": [2, 2, 3, 4, 4, 4, 6]
    }
    """
    all_distances = {}
    for pair in pairz:
        sel1 = getCenterOfMassCoordinates(universe, f"protein and resid {pair[0]}")
        sel2 = getCenterOfMassCoordinates(universe, f"protein and resid {pair[1]}")
        groupLabel = f"{pair[0]}-{pair[1]}"
        p_distances = getPairWiseDistance(sel1, sel2)
        all_distances[groupLabel] = p_distances

    return all_distances


def formatForPlotting(moved_pdistances, stat_pdistances, neu_pdistances):

    df1 = pd.DataFrame(moved_pdistances)
    df1 = df1.melt(var_name="Residue Pair", value_name="Distance")
    df1["State"] = "Transitioned"

    df2 = pd.DataFrame(stat_pdistances)
    df2 = df2.melt(var_name="Residue Pair", value_name="Distance")
    df2["State"] = "Stationary"

    df3 = pd.DataFrame(neu_pdistances)
    df3 = df3.melt(var_name="Residue Pair", value_name="Distance")
    df3["State"] = "Neutral"

    full_df = pd.concat([df1, df2, df3], ignore_index=True)

    return full_df


def plotSimplePairwiseDistances(trad, stad, nmd):
    rc('font', **{'family': 'serif', 'serif': ['Arial']})
    rc('text', usetex = True)
    fig, ax = plt.subplots(figsize = (15, 8))

    residue_pair = "163-215"

    sns.kdeplot(x=trad[residue_pair], color="#d0bbff", linewidth = 3)
    sns.kdeplot(x=stad[residue_pair], color = "#b9f2f0", linewidth = 3)
    sns.kdeplot(x=nmd[residue_pair], color="#ff9f9b", linewidth = 3)

    figure = ax.get_figure()
    figure.savefig('pairwises_r.pdf', dpi = 500)
    ax.set(title = 'Pairwise Residue Fluctuations')
    plt.show()

if __name__ == '__main__':
    matplotlib.rcParams['font.family'] = 'Arial'
    plt.style.use(['science', 'ieee', 'no-latex'])

    psf1 = '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/stripped.psf'
    dcd1 = '/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_and_aligned/smallChargedMembrane-ions-wats-00004-stripped.dcd'

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
    datafr = formatForPlotting(tra_pdist, stat_pdist, nm_pdist)
    plotSimplePairwiseDistances(tra_pdist, stat_pdist, nm_pdist)


