from typing import Any
from collections import Counter

class MDStats:
    """A collection of functions designed to calculate several different statistics involved in Molecular Dynamics Simulations."""
    def __init__(self):
        pass

    @staticmethod
    def countDistributions(collection: list, tolerance: int) -> list:
        """Takes in a list of same shape collections of numerical parameters and returns a list of tuples in the form: (VALUES_OF_PARAMETER, COUNTS)"""
        #[(x, y, z), (x, y, z)] -> 3 distributions for x, y, and z
        num_parameters = len(collection[0])
        distributions = []

        for ind in range(num_parameters):
            #Get parameter values
            parameter_vals = []
            counts = {}
            sorted_collection = sorted(collection, key=lambda group: group[ind])
            for data_group in sorted_collection:
                parameter_vals.append(data_group[ind])

            #Get bins
            binned_parameter_vals = [(round(x, tolerance)) for x in parameter_vals]

            #Calculate distribution
            for binned_parameter in binned_parameter_vals:
                counts[binned_parameter] = counts.get(binned_parameter, 0) + 1

            distributions.append(counts)
        return distributions


    @staticmethod
    def countSetDistribution(collection: list, tolerance: int):
        counts = {}

        # Bin parameter sets
        parameters = [(round(euler_set[0], tolerance), round(euler_set[1], tolerance), round(euler_set[2], tolerance)) for euler_set in collection if euler_set != (0.0, -0.0, 0.0)]

        # Calculate Distribution
        for binned_parameter in parameters:
            counts[binned_parameter] = counts.get(binned_parameter, 0) + 1

        return counts


    @staticmethod
    def findModes(collection: list, num_modes: int) -> tuple:
        """Takes in a distribution and finds the modes of the dataset."""
        counts = Counter(collection)
        multimodes = [x[0] for x in counts.most_common(num_modes)]

        return multimodes


    @staticmethod
    def removeDuplicates(repetitive_items: Any) -> list:
        """Removes duplicates yet maintains order"""
        seen = {}
        duplicates_removed = []
        for item in repetitive_items:
            if item not in seen:
                duplicates_removed.append(item)
            seen[item] = None

        return duplicates_removed

