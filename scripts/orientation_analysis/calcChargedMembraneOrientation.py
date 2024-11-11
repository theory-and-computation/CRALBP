import sys
sys.path.append('/Users/danielsantos/Downloads/CRALBP_analysis/modules')

import numpy as np
import pandas as pd
import MDdistance as dist
import MDStats as mdstats
import Parse as parse
import TCLcalls as tcl
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib.ticker as mtick


def calcRotMatrix(psf: str, base_dcd: str, outfile: str, refpdb: str, traj_num: int) -> None:
    """Calls tcl script that writes out a rotation matrix describing the protein's orientation in each frame."""

    for i in range(traj_num):
        if i + 4 < 10:
            dcd = base_dcd + f"-0000{i+4}-stripped.dcd"
        else:
            dcd = base_dcd + f"-000{i+4}-stripped.dcd"

        tcl.TCLcalls.executeGetRotationalMatrices(psf, dcd, outfile, refpdb)


def parseRotMatrix(filename) -> list:
    """Parses rotation matrix temporary text file and transforms each into a numpy array"""
    filelines = parse.Parse.readFile(filename)
    rotation_matrices = parse.Parse.parseRotMatrix(filelines)

    return rotation_matrices


def calcEulerAngles(rot_matrices) -> list:
    """Given a list of 3 x 3 rotation matrices, extracts all yaw, pitch, and roll angles."""
    euler_angles = []
    for rot_matrix in rot_matrices:
        euler_angles.append(dist.MDdistance.getEulerAnglesFromRotMatrix(rot_matrix))
    return euler_angles


def calcDistributions(euler_angles: tuple) -> pd.DataFrame:
    """Calculates distribution of the euler angles throughout the trajectory and formats distribution data into a pandas Dataframe for plotting with Seaborn"""
    distributions = mdstats.MDStats.countDistributions(euler_angles, 1)
    yaw_dat = pd.DataFrame(distributions[0], index=range(len(distributions[0])))
    pitch_dat = pd.DataFrame(distributions[1], index=range(len(distributions[1])))
    roll_dat = pd.DataFrame(distributions[2], index=range(len(distributions[2])))

    return yaw_dat, pitch_dat, roll_dat


def plotDistributions(dat, angle: str):

    # Interpolate the radii
    theta = list(dat.columns.values)
    radii = list(dat.iloc[1])
    interp_func = interpolate.interp1d(theta, radii, fill_value = 'extrapolate', kind = 'cubic')
    new_radii = interp_func(theta)

    # Get probabilities
    width = 2 * np.pi / len(theta)
    prob = np.array(new_radii) / new_radii.sum()

    #Set-up:
    fig = plt.figure()
    ax = fig.add_subplot(111, polar = True)

    if angle == 'yaw':

        #Figure settings
        bars = ax.bar(theta, prob, width = width, bottom = 0.0, color = "#ff9f9b")
        _max = 0.045
        ax.set_rorigin(-(_max / 6))
        ax.set_rmax(_max)
        ax.set_rticks([0.015, 0.030])
        ax.set_rlabel_position(0)
        ax.tick_params(axis='y', labelsize = 12)
        ax.tick_params(axis = 'x', labelsize = 15, pad=10)
        ax.grid(True)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals = 0))
        ax.tick_params(axis = 'y', labelsize = 12)
        ax.tick_params(axis = 'x', labelsize = 15, pad = 10)
        ax.set_title('Alpha Angle Distributions', va = 'bottom')

        #Save
        figure = ax.get_figure()
        figure.savefig('/Users/danielsantos/Downloads/figures/orientation_CRALBP/yaw_distribution.pdf', bbox_inches='tight', dpi = 500)

    elif angle == 'pitch':

        #Figure settings
        bars = ax.bar(theta, prob, width = width, bottom = 0.0, color = "#d0bbff")
        _max = 0.18
        ax.set_rorigin(-(_max / 6))
        ax.set_rmax(_max)
        ax.set_rticks([0.06, 0.12])
        ax.set_rlabel_position(0)
        ax.tick_params(axis = 'y', labelsize = 12)
        ax.tick_params(axis = 'x', labelsize = 15, pad=10)
        ax.grid(True)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals = 1))
        ax.set_title('Beta Angle Distributions', va = 'bottom')

        #Save
        figure = ax.get_figure()
        figure.savefig('/Users/danielsantos/Downloads/figures/orientation_CRALBP/pitch_distribution.pdf',bbox_inches='tight', dpi = 500)

    elif angle == 'roll':

        #Figure settings
        bars = ax.bar(theta, prob, width = width, bottom = 0.0, color = "#b9f2f0")
        _max = 0.16
        ax.set_rorigin(-(_max / 6))
        ax.set_rmax(_max)
        ax.set_rticks([0.05, 0.10])
        ax.set_rlabel_position(0)
        ax.tick_params(axis = 'y', labelsize = 12)
        ax.tick_params(axis = 'x', labelsize = 15, pad = 10)
        ax.grid(True)
        ax.yaxis.set_major_formatter(mtick.PercentFormatter(1, decimals = 0))
        ax.set_title('Gamma Angle Distributions', va = 'bottom')

        #Save
        figure = ax.get_figure()
        figure.savefig('/Users/danielsantos/Downloads/figures/orientation_CRALBP/roll_distribution.pdf', bbox_inches = 'tight', dpi = 500)

    plt.show()


def visualizePrefOrientation(inpdb, euler_angle_set, outfile_set) -> None:
    """Calls tcl script that moves the protein according to its preferred position and writes to a pdb."""

    if len(euler_angle_set) != len(outfile_set):
        raise Exception("Number of outfiles must match number of preferred orientations")

    for ind in range(len(euler_angle_set)):
        yaw = np.degrees(euler_angle_set[ind][0])
        pitch = np.degrees(euler_angle_set[ind][1])
        roll = np.degrees(euler_angle_set[ind][2])

        tcl.TCLcalls.executeRotateProteinToPreferred(inpdb, yaw, pitch, roll, outfile_set[ind])



if __name__ == '__main__':


    PSF = "/Users/danielsantos/Downloads/CRALBP_analysis/SMC_STRIPPED_2/stripped.psf"
    BASE_DCD = "/Users/danielsantos/Downloads/CRALBP_analysis/SMC_STRIPPED_2/smallChargedMembrane_2-ions-wats"
    REFPDB = "/Users/danielsantos/Downloads/CRALBP_analysis/SMC_STRIPPED_2/first_prot.oriented.pdb"

    NUMDCD = 30

    OUTTXT = "/Users/danielsantos/Downloads/CRALBP_analysis/SMC_stripped_2/rotMat.txt"
    OUTPDBS = ["/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode1.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode2.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode3.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode4.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode5.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode6.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode7.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode8.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode9.pdb",
              "/Users/danielsantos/Downloads/CRALBP_analysis/misc/SMC_2_mode10.pdb"]


    calcRotMatrix(PSF, BASE_DCD, OUTTXT, REFPDB, NUMDCD)
    rots = parseRotMatrix(OUTTXT)
    euler_sets = calcEulerAngles(rots)
    set_distribution = mdstats.MDStats.countSetDistribution(euler_sets, 0)

    angles = ['yaw', 'pitch', 'roll']
    dataframes = list(calcDistributions(euler_sets))

    for ind in range(len(angles)):
        plotDistributions(dataframes[ind], angles[ind])

    mods = mdstats.MDStats.findModes(set_distribution, 10)
    visualizePrefOrientation(REFPDB, mods[:], OUTPDBS)

