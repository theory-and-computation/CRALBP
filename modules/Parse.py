import numpy as np

class Parse:
    """A collection of functions designed to parse text files."""
    def __init__(self):
        pass

    @staticmethod
    def readFile(filename) -> list:
        """Reads from a text file."""
        with open(filename, 'r') as fileobj:
            lines = fileobj.readlines()
            return lines

    @staticmethod
    def parseRotMatrix(filelines) -> list:
        """Parses 3x3 rotation matrices from text file containing 4x4 VMD transformation matrices"""
        rotation_matrices = []
        for line in filelines:
            no_ends = line[1:-2]
            split_mat = no_ends.split("} {")
            file_matrix = []
            for row in split_mat:
                split_row = row.split(" ")
                file_matrix.append(split_row)

            rot_mat = np.array([[file_matrix[0][0], file_matrix[0][1], file_matrix[0][2]],
                                [file_matrix[1][0], file_matrix[1][1], file_matrix[1][2]],
                                [file_matrix[2][0], file_matrix[2][1], file_matrix[2][2]]])

            rotation_matrices.append(rot_mat)

        return rotation_matrices


