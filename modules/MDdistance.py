import MDAnalysis as mda
import numpy as np
import math

class MDdistance:
    """A collection of functions useful for calculating certain parameters in a 3D Euclidean space"""
    def __init__(self):
        pass

    @staticmethod
    def getEulerAnglesFromRotMatrix(rot_matrix) -> (float, float, float):
        """Given a 3 x 3 rotation matrix, the Euler angles: yaw, pitch, and roll are calculated. Note: RETURNS ANGLES IN RADIANS"""
        pos00, pos01, pos02 = float(rot_matrix[0, 0]), float(rot_matrix[0, 1]), float(rot_matrix[0, 2])
        pos10, pos11, pos12 = float(rot_matrix[1, 0]), float(rot_matrix[1, 1]), float(rot_matrix[1, 2])
        pos20, pos21, pos22 = float(rot_matrix[2, 0]), float(rot_matrix[2, 1]), float(rot_matrix[2, 2])

        yaw = np.mod(np.arctan2(pos10, pos00), 2 * np.pi)
        pitch = np.mod(np.arctan2(-pos20, np.sqrt(pos21**2 + pos22**2)), 2 * np.pi)
        roll = np.mod(np.arctan2(pos21, pos22), 2 * np.pi)

        return yaw, pitch, roll


    @staticmethod
    def getRotMatrixFromEulerAngles(yaw, pitch, roll) -> np.array:
        """Given the Euler angles: yaw, pitch, and roll, a 3 x 3 rotation matrix is calculated."""
        yaw = np.radians(yaw)
        pitch = np.radians(pitch)
        roll = np.radians(roll)

        Rx = np.array([[1, 0, 0], [0, np.cos(roll), -np.sin(roll)], [0, np.sin(roll), np.cos(roll)]])
        Ry = np.array([[np.cos(pitch), 0, np.sin(pitch)], [0, 1, 0], [-np.sin(pitch), 0, np.cos(pitch)]])
        Rz = np.array([[np.cos(yaw), -np.sin(yaw), 0], [np.sin(yaw), np.cos(yaw), 0], [0, 0, 1]])

        R = np.dot(np.dot(Rx, Ry), Rz)

        return R


    @staticmethod
    def getCOMCoordinates(*, universe: mda.Universe = None, selection_command: str = "", frame_count: int = 0) -> dict[int, str]:
        """Obtain the coordinates for a given selection using MD analysis."""
        selection_coordinates = {}

        if selection_command != "" and frame_count != 0:
            try:
                select = universe.select_atoms(selection_command)
                for frame_index in range(frame_count):

                    universe.trajectory[frame_index]
                    selection_COM = select.center_of_mass()
                    selection_coordinates[frame_index] = selection_COM

                return selection_coordinates
            except:
                print("Could not make selection.")

        return {}


    @staticmethod
    def calculate3DDist(*, ThreeDCoords1: tuple, ThreeDCoords2: tuple, box_dimensions: tuple, adjustForMemb=False):
        """Takes in 2 three-dimensional coordinate sets and calculates the distance from coordinate 1 to coordinate 2."""
        x_diff = abs(ThreeDCoords2[0] - ThreeDCoords1[0])
        y_diff = abs(ThreeDCoords2[1] - ThreeDCoords1[1])
        z_diff = abs(ThreeDCoords2[2] - ThreeDCoords1[2])

        if adjustForMemb:
            diffs = MDdistance.adjustForPBC(x_diff, y_diff, z_diff, box_dimensions)
        else:
            diffs = (x_diff, y_diff, z_diff)

        return math.sqrt((diffs[0]**2) + (diffs[1]**2) + (diffs[2]**2))


    @staticmethod
    def adjustForPBC(x_diff: float, y_diff: float, z_diff: float, box_dimensions: tuple) -> tuple:
        """Adjusts coordinate differences according to the periodic boundary. Accounts for jumps."""
        cutoff_x = box_dimensions[0] / 2
        cutoff_y = box_dimensions[1] / 2
        cutoff_z = box_dimensions[2] / 2

        new_x_diff = x_diff
        new_y_diff = y_diff
        new_z_diff = z_diff

        if x_diff > cutoff_x:
            new_x_diff = box_dimensions[0] - x_diff
        if y_diff > cutoff_y:
            new_y_diff = box_dimensions[1] - y_diff
        if z_diff > cutoff_z:
            new_z_diff = box_dimensions[2] - z_diff

        return new_x_diff, new_y_diff, new_z_diff


    @staticmethod
    def adjustForMembraneDiffusion(threeDCoordSet: tuple , membCoordSet) -> tuple:
        """Takes in an object three-dimensional coordinate set and a membrane center of mass coordinte set and
        return the object's new coordinate values with its z now relative to the membrane's center of mass."""
        new_coord_set = (threeDCoordSet[0], threeDCoordSet[1], threeDCoordSet[2] - membCoordSet[2])
        return new_coord_set
