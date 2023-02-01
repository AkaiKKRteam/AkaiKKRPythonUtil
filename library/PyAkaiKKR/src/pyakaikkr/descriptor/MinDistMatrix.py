# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import numpy as np
from pymatgen.analysis.structure_analyzer import VoronoiConnectivity
from pymatgen.core import Structure

def min_dist_matrix(structure: Structure) -> float:
    """make minimum distance of the structure from VoronoiConnectivity(structure)

    Args:
        structure (pymatgen.core.structure): structure

    Returns:
        float: minium distance of the distance matrix
    """
    n = len(structure.sites)
    min_dist_matrix = np.zeros((n, n))
    vorcon = VoronoiConnectivity(structure)
    for conn in vorcon.get_connections():
        i = conn[0]
        j = conn[1]
        d = conn[2]
        min_dist_matrix[i, j] = d
    # print(min_dist_matrix)
    return min_dist_matrix

