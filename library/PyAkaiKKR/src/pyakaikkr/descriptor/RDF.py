# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.



import shutil
from pymatgen.core.sites import PeriodicSite
from pymatgen.core import Structure

from pymatgen.analysis.structure_matcher import StructureMatcher
import json
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifParser

import sys
import subprocess
import hashlib
import os

import numpy as np
import matplotlib.pyplot as plt



sys.path.append("/home/kino/tmp/AkaiKKRPythonUtil/util")
if "test_script" not in sys.modules:
    from pyakaikkr import AkaikkrJob
    from pyakaikkr.Error import *
    from pyakaikkr import StructureSpeciesConverter


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


class RDFConverter:
    def __init__(self):
        pass

    def make_nndistance(self, structure: Structure, rcut=4.0,) -> dict:
        """make nndistance

        Args:
            structure (pymatgen.core.structure): struture
            rcut (float, optional): cut off distance. Defaults to 4.0.
        Returns:
            dict: 'sum': all the nndistance, 'locals': nndistance set for each site, 'atoms': each site
        """
        atoms = np.array([site.species_string for site in structure])

        local_nndistance = []
        for site in structure.sites:
            neighbor = structure.get_neighbors(site, r=rcut)
            nndistance = []
            for nn in neighbor:
                nndistance.append(nn[1])
            local_nndistance.append(nndistance)
        nndistance = []
        for d in local_nndistance:
            nndistance.extend(d)

        return {'all': nndistance, 'locals': local_nndistance, 'atoms': atoms}

    def transform(self, structure: Structure = None, nndistance=None,
                  distancetype='all',
                  range=None, rcut=4.0, bins=200):
        """transform structure to RDF

        structure or nndistance is necessary.

        If nndistance is None and structure is given, RDF is made of 

        Args:
            structure (Structure, optional): structure 
            nndistance (dict, optional): nndistance. Defaults to None.
            distancetype (str, optional): distance type. Defualts to 'all'.
            range (int, optional): min and man vlaues of bins. Defalts to None
            rcut (float, optional): distance cutoff value. Defaults to 4.0.
            bins (int, optional): number of bins in the histogram. Defaults to 200.

        Returns:
            np.ndarray, np.ndarray:  hist, edges
        """
        if structure is not None and nndistance is None:
            nndistance = self.make_nndistance(structure, rcut=rcut)
        nn = nndistance[distancetype]

        if distancetype == 'all':
            if range is None:
                minvalue = np.min(nn)
                maxvalue = np.max(nn)
                range = (minvalue, maxvalue)
            hist, edges = np.histogram(nn, bins=bins, range=range)
            return hist, edges
        elif distancetype == 'locals':

            raise AttributeError(
                'distancetype={} not supported now'.format(distancetype))


def RDFHelper_fig(hist, edges, ax=None, directory="output"):
    """plot input and output RDF

    Args:
        hist (nd.array): RDF histogram
        edge (nd.array): RDF edges
        directory (str): directory to save png
        bins (int, optional): number of bins. Defaults to 200.

    """

    binvalue = (edges[:-1] + edges[1:])*0.5
    width = binvalue[1]-binvalue[0]

    ax_org = ax
    if ax_org is None:
        fig, ax = plt.subplots()
    ax.bar(binvalue, hist, width=width)

    if ax_org is None and directory is not None:
        os.makedirs(directory, exist_ok=True)
        fig.tight_layout()
        filename = os.path.join(directory, "rdf.png")
        fig.savefig(filename)
        print("saved to", filename)
        fig.clf()
        plt.close(fig)


if __name__ == "__main__":
    def main():

        prefix = "/home/kino/tmp/AkaiKKRPythonUtil/util/cif2kkr_test_script"
        ciffile = "MaterialsLibrary/AtomWorkData/spcno_series/spcno_088_AgO.cif"
        ciffilepath = os.path.join(prefix, ciffile)

        structure = Structure.from_file(ciffilepath)
        rdfconverter = RDFConverter()

        hist, edge = rdfconverter.transform(structure)

        RDFHelper_fig(hist, edge)

    main()
