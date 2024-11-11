import sys
sys.path.append('/Users/danielsantos/Downloads/CRALBP_analysis/modules')

import MDAnalysis as mda
import math
import numpy as npy
import MDdistance as dist

#Three-Dimensional.

def initializeUniverse(DCD, PSF) -> mda.Universe:
    """Initializes the universe object."""
    uni = mda.Universe(PSF, DCD)

    return uni


def getCenterOfMassCoordinates(universe) -> (dict, int):
    """Calculates center of mass coordinates for a protein in a trajectory file."""
    total_frames = len(universe.trajectory)
    protein = universe.select_atoms('protein')

    protein_COM_coordinates = {}

    for frame_index in range(total_frames):

        universe.trajectory[frame_index]
        protein_COM = protein.center_of_mass()

        protein_COM_coordinates[frame_index] = protein_COM

    return protein_COM_coordinates, total_frames


def getMembraneCenterOfMassCoordinates(universe) -> (dict, int):

    total_frames = len(universe.trajectory)
    membrane = universe.select_atoms('segid MEMB')

    membrane_COM_coordinates = {}

    for frame_index in range(total_frames):
        universe.trajectory[frame_index]
        membrane_COM = membrane.center_of_mass()

        membrane_COM_coordinates[frame_index] = membrane_COM

    return membrane_COM_coordinates, total_frames


def runMSD(COM_coordinates, MEMB_coordinates, total_frame_count, box_size) -> (list, list):
    """Runs Mean Square Deviation calculations."""
    lags = []
    msds = []
    standard_devs = []
    standard_errors = []

    for lag_time in range(10000):
        print(f'lag time: {lag_time}')
        squared_displacements = []

        for frame_index in range(total_frame_count):

            if frame_index + lag_time < total_frame_count:

                protein_COM_initial = COM_coordinates[frame_index]
                protein_COM_final = COM_coordinates[frame_index + lag_time]
                MEMB_COM_initial = MEMB_coordinates[frame_index]
                MEMB_COM_final = MEMB_coordinates[frame_index + lag_time]
                box_dims = box_size[frame_index + lag_time]

                adjusted_protCOM_initial = dist.MDdistance.adjustForMembraneDiffusion(protein_COM_initial, MEMB_COM_initial)
                adjusted_protCOM_final = dist.MDdistance.adjustForMembraneDiffusion(protein_COM_final, MEMB_COM_final)

                displacement = dist.MDdistance.calculate3DDist(ThreeDCoords1=adjusted_protCOM_initial, ThreeDCoords2=adjusted_protCOM_final, box_dimensions=box_dims, adjustForMemb=True)

                squared_displacement = displacement ** 2
                squared_displacements.append(squared_displacement)


        #Calculating standard deviation, standard error, and storing data.
        frame_MSD = sum(squared_displacements) / len(squared_displacements)

        std_dev = calcStandardDeviation(squared_displacements)
        std_error = calcStandardError(squared_displacements, std_dev)

        lags.append(lag_time)
        msds.append(frame_MSD)
        standard_devs.append(std_dev)
        standard_errors.append(std_error)

    return lags, msds, standard_devs, standard_errors


def calcSquaredDisplacement(x_diff, y_diff, z_diff):

    _squared_x_difference = math.pow(x_diff, 2)
    _squared_y_difference = math.pow(y_diff, 2)
    _squared_z_difference = math.pow(z_diff, 2)

    _displacement = math.sqrt(_squared_x_difference + _squared_y_difference + _squared_z_difference)
    _squared_displacement = math.pow(_displacement, 2)

    return _squared_displacement


def calcStandardDeviation(dataset: list):

    if len(dataset) < 2:

        standard_deviation = 0
        return standard_deviation

    else:

        dataset_mean = sum(dataset) / len(dataset)
        squared_diff = [(val - dataset_mean) ** 2 for val in dataset]
        variance = sum(squared_diff) / (len(dataset) - 1)

        standard_deviation = variance ** 0.5

        return standard_deviation


def calcStandardError(dataset: list, standard_dev: float):

    if standard_dev == 0:
        standard_error = 0

    else:
        standard_error = standard_dev / math.sqrt(len(dataset))

    return standard_error


def checkPBC(x_coordinate_info: tuple, y_coordinate_info: tuple, z_coordinate_info: tuple, box_dimensions: tuple):

    box_length = box_dimensions[0]
    box_width = box_dimensions[1]
    box_height = box_dimensions[2]

    x_diff = x_coordinate_info[2]
    y_diff = y_coordinate_info[2]
    z_diff = z_coordinate_info[2]

    cutoff_x = box_length / 2
    cutoff_y = box_width / 2
    cutoff_z = box_height / 2

    new_x_diff = x_diff
    new_y_diff = y_diff
    new_z_diff = z_diff


    if x_diff > cutoff_x:
        new_x_diff = box_length - x_diff
    if y_diff > cutoff_y:
        new_y_diff = box_width - y_diff
    if z_diff > cutoff_z:
        new_z_diff = box_height - z_diff

    return new_x_diff, new_y_diff, new_z_diff


def saveData(charged_membrane_data, neutral_membrane_data, CM_stdev, NM_stdev, CM_std_error, NM_std_error):

    both_datasets = [charged_membrane_data, neutral_membrane_data, CM_stdev, NM_stdev, CM_std_error, NM_std_error]
    numpified = npy.array(both_datasets)

    npy.save('/Users/danielsantos/Downloads/CRALBP/save_files/MSD.npy', numpified)


def readBoxFile(b_file):

    with open(b_file, 'r') as fileobj:
        lines = fileobj.readlines()
    return lines


def parseBoxLines(file_lines):

    box_vals = []

    for line in file_lines:

        stripped = line.strip()
        removed_braces_open = stripped.replace('{', '')
        removed_braces = removed_braces_open.replace('}', '')

        dimensions = removed_braces.split()

        box_length = float(dimensions[0])
        box_width = float(dimensions[1])
        box_height = float(dimensions[2])

        _dims = (box_length, box_width, box_height)
        box_vals.append(_dims)

    return box_vals


if __name__ == '__main__':

    dcd = '../RSMC_stripped/rsmc.dcd'
    psf = '../RSMC_stripped/stripped.psf'
    CM_BOX_FILE = './RSMC_PBC_SIZE.txt'

    dcd2 = '../all_stripped/FNM.dcd'
    psf2 = '../all_stripped/stripped.psf'
    NM_BOX_FILE = './FNM_PBC_SIZE.txt'

    u = initializeUniverse(dcd, psf)
    second_u = initializeUniverse(dcd2, psf2)

    CM_membrane_COM_coords_uni_one, frame1 = getMembraneCenterOfMassCoordinates(u)
    NM_membrane_COM_coords_uni_two, frame2 = getMembraneCenterOfMassCoordinates(second_u)

    COM_coords_uni_one, frame_count_uni_one = getCenterOfMassCoordinates(u)
    COM_coords_uni_two, frame_count_uni_two = getCenterOfMassCoordinates(second_u)

    CM_BOX_FILE_LINES = readBoxFile(CM_BOX_FILE)
    CM_BOX = parseBoxLines(CM_BOX_FILE_LINES)

    NM_BOX_FILE_LINES = readBoxFile(NM_BOX_FILE)
    NM_BOX = parseBoxLines(NM_BOX_FILE_LINES)

    lag, temp_msd, stdevs, stderror = runMSD(COM_coords_uni_one, CM_membrane_COM_coords_uni_one, frame_count_uni_one, CM_BOX)
    lag2, temp_msd2, stdevs2, stderror2 = runMSD(COM_coords_uni_two, NM_membrane_COM_coords_uni_two, frame_count_uni_two, NM_BOX)

    saveData(temp_msd, temp_msd2, stdevs, stdevs2, stderror, stderror2)

