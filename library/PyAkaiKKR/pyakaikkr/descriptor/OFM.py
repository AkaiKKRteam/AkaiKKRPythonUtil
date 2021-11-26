# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import re
import numpy as np
import os
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.core import Structure
from pymatgen.core.periodic_table import Element
import pandas as pd
from copy import deepcopy
import seaborn as sns
import matplotlib.pyplot as plt
from typing import List

class OFMElementRepresentation:
    @staticmethod
    def as_electronic_structure(name=None,
                                show_detail=False,
                                **kwarg) -> np.ndarray:
        """
        generate one-hot representation for a element, e.g, si = [0.0, 1.0, 0.0, 0.0, ...]
        :param name: element symbol
        """
        if name is None:
            name = "H"

        if "is_including_row" in kwarg:
            is_including_row = kwarg["is_including_row"]
        else:
            is_including_row = False

        element = Element(name)

        general_electron_subshells = ['s1', 's2', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6',
                                      'd1', 'd2', 'd3', 'd4', 'd5', 'd6', 'd7', 'd8', 'd9', 'd10',
                                      'f1', 'f2', 'f3', 'f4', 'f5', 'f6', 'f7',
                                      'f8', 'f9', 'f10', 'f11', 'f12', 'f13', 'f14']

        general_element_electronic = {}
        for key in general_electron_subshells:
            general_element_electronic[key] = 0.0

        if name == 'H':
            element_electronic_structure = ['s1']
        elif name == 'He':
            element_electronic_structure = ['s2']
        else:
            element_electronic_structure = [''.join(pair) for pair in re.findall("\.\d(\w+)(\d+)",
                                                                                 element.electronic_structure)]
        if show_detail:
            print(element_electronic_structure)
        for eletron_subshell in element_electronic_structure:
            general_element_electronic[eletron_subshell] = 1.0

        localofm = np.array([general_element_electronic[key]
                            for key in general_electron_subshells])

        if is_including_row:
            for irow in range(1, 10):
                general_electron_subshells.append("row{}".format(irow))
            localofm_row = np.zeros(9)
            irow = element.row-1
            localofm_row[irow] = 1.0
            if show_detail:
                print(localofm_row)
            localofm = np.hstack([localofm, localofm_row])

        return localofm, general_electron_subshells

    @staticmethod
    def as_periodic_table(name=None,
                          show_detail=False, **kwarg) -> np.ndarray:
        """
        generate one-hot representation for a element, e.g, si = [0.0, 1.0, 0.0, 0.0, ...]
        :param name: element symbol
        """
        if name is None:
            name = "H"

        if "is_including_row" in kwarg:
            is_including_row = kwarg["is_including_row"]
        else:
            is_including_row = False

        element = Element(name)

        localofm = np.zeros(18)
        columns = []
        for i in range(1, 19):
            columns.append("group{}".format(i))

        localofm[element.group-1] = 1.0
        if show_detail:
            print(element.group, element.row)
            print(localofm)

        if is_including_row:
            for irow in range(1, 10):
                columns.append("row{}".format(irow))
            localofm_row = np.zeros(9)
            irow = element.row-1
            localofm_row[irow] = 1.0
            if show_detail:
                print(localofm_row)
            localofm = np.hstack([localofm, localofm_row])

        return localofm, columns

    @staticmethod
    def as_name(name=None,
                show_detail=False, **kwarg) -> np.ndarray:
        """
        generate one-hot representation for a element, e.g, si = [0.0, 1.0, 0.0, 0.0, ...]
        :param name: element symbol
        """
        names = kwarg["names"]

        if name is None:
            name = names[0]

        localofm = np.zeros(len(names))
        i = names.index(name)
        if show_detail:
            print(name, i, names)
        localofm[i] = 1.0

        return localofm, names


class OFMConverter:
    def __init__(self, is_ofm1=True, is_including_d=True,
                 representation="electronic_structure",
                 **kwarg):
        """
        If element_representation=="electronic_structure" or element_representation=="periodic_table",
            kwarg["is_including_row"] can define whether including the row number of the periodic table.
        If element_representation=="name",
            kwarg["name"] must be defined. It is a list of all names.


        Args:
            is_ofm1 (bool, optional): OFM1. Default to True.
            is_including_d (bool, optional): including 1/d contribution. Default to True.
            representation (str, optional): element_representation. Default to "electronic_structure".

        """
        self.is_ofm1 = is_ofm1
        self.is_including_d = is_including_d
        self.representation = representation

        if representation == "electronic_structure":
            self.element_representation = OFMElementRepresentation.as_electronic_structure
        elif representation == "periodic_table":
            self.element_representation = OFMElementRepresentation.as_periodic_table
        elif representation == "name":
            self.element_representation = OFMElementRepresentation.as_name
            # arg["name"] is necessary
            if "names" not in kwarg:
                raise AttributeError(
                    "'names' values are necessary in kwarg in representation {}".format(representation))
        else:
            raise AttributeError(
                "unknown representation {}".format(representation))

        self.kwarg = kwarg

    def transform(self, struct: Structure) -> dict:
        """
        Generate OFM descriptor
        :param struct: pymatgen Structure object
        """

        is_ofm1 = self.is_ofm1
        is_including_d = self.is_including_d
        element_representation = self.element_representation
        kwarg = self.kwarg

        atoms = np.array([site.species_string for site in struct])

        vector, index = element_representation(atom=None, **kwarg)

        if is_ofm1:
            columns = deepcopy(index)
            columns.insert(0, "1")

        local_xyz = []
        local_orbital_field_matrices = []
        for i_atom, atom in enumerate(atoms):

            coordinator_finder = VoronoiNN(cutoff=10.0)
            neighbors = coordinator_finder.get_nn_info(
                structure=struct, n=i_atom)

            site = struct[i_atom]
            center_vector, _ = element_representation(atom, **kwarg)
            env_vector = np.zeros_like(vector)

            atom_xyz = [atom]
            coords_xyz = [site.coords]

            for nn in neighbors:
                site_x = nn['site']
                w = nn['weight']
                site_x_label = site_x.species_string
                atom_xyz += [site_x_label]
                coords_xyz += [site_x.coords]
                neigh_vector, _ = element_representation(site_x_label, **kwarg)
                d = np.sqrt(np.sum((site.coords - site_x.coords)**2))
                if is_including_d:
                    env_vector += neigh_vector * w / d
                else:
                    env_vector += neigh_vector * w

            if is_ofm1:
                env_vector = np.concatenate(([1.0], env_vector))

            local_matrix = center_vector[None, :] * env_vector[:, None]
            if False:
                print(center_vector)
                print(env_vector)
                print(local_matrix)
                raise

            # to 1D
            local_matrix = np.ravel(local_matrix)
            local_orbital_field_matrices.append(local_matrix)
            local_xyz.append({"atoms": np.array(atom_xyz),
                             "coords": np.array(coords_xyz)})

        local_orbital_field_matrices = np.array(local_orbital_field_matrices)
        material_descriptor = np.mean(local_orbital_field_matrices, axis=0)

        return {'mean': material_descriptor,
                'locals': local_orbital_field_matrices,
                'atoms': atoms,
                'columns': columns,
                'index': index,
                "local_xyz": local_xyz}


def OFMHelper_structure2names(struct: Structure) -> List[str]:
    names = []
    for site in struct.sites:
        names.append(str(site.species_string))
    names.sort()
    names = list(set(names))
    return names


def OFMHelper_local_df(ofm: dict, i) -> pd.DataFrame:
    index = ofm["index"]
    columns = ofm["columns"]
    n1 = len(index)
    n2 = len(columns)
    localofm = ofm["locals"][i]
    localofm = localofm.reshape((n2, n1)).T
    df = pd.DataFrame(localofm, index=index, columns=columns)
    return df


def OFMHelper_mean_df(ofm: dict) -> pd.DataFrame:
    index = ofm["index"]
    columns = ofm["columns"]
    n1 = len(index)
    n2 = len(columns)
    localofm = ofm["mean"]
    localofm = localofm.reshape((n2, n1)).T
    df = pd.DataFrame(localofm, index=index, columns=columns)
    return df


def OFMHelper_fig(ofm: dict, action=["mean", "locals"], annot=False):
    """generate OFM plots

    If "mean" in action, plot ofm["mean"].
    If "locals" in action, plot ofm["locals].

    Args:
        ofm (dict): output of OFMCovnerter
        action ([str], optional): actions to do. Defaults to ["mean", "locals"].
    """
    outputdir = "output"

    if "mean" in action:
        os.makedirs(outputdir, exist_ok=True)

        df = OFMHelper_mean_df(ofm)
        fig, ax = plt.subplots()
        sns.heatmap(df,  annot=annot, ax=ax)
        pngfile = "mean.png"
        pngfilepath = os.path.join(outputdir, pngfile)
        fig.tight_layout()
        plt.savefig(pngfilepath)
        fig.clf()
        plt.close(fig)
        print(pngfilepath, "is made.")

    if "locals" in action:
        os.makedirs(outputdir, exist_ok=True)

        for i, name in enumerate(ofm["atoms"]):
            df = OFMHelper_local_df(ofm, i)
            fig, ax = plt.subplots()
            ax.set_title("center={}".format(name))
            sns.heatmap(df, annot=annot, ax=ax)
            pngfile = "local{}.png".format(i)
            pngfilepath = os.path.join(outputdir, pngfile)
            fig.tight_layout()
            plt.savefig(pngfilepath)
            fig.clf()
            plt.close(fig)
            print(pngfilepath, "is made.")


if __name__ == "__main__":

    def get_ciffiles(prefix="/home/kino/kino/kit/AkaiKKRPythonUtil/util/cif2kkr_test_script"):
        if True:
            files = [
                "MaterialsLibrary/AtomWorkData/spcno_series/spcno_042_Bi5FeO15Ti3.cif"]
        if True:
            files = [
                "MaterialsLibrary/Intermetallics/Eu3Ni4Ga4 (I-43m).cif"]
        filepaths = []
        for file in files:
            filepaths.append(os.path.join(prefix, file))
        return filepaths

    def main():

        filepaths = get_ciffiles()

        file = filepaths[0]
        struct = Structure.from_file(file)

        representation = "electronic_structure"
        # representation = "periodic_table"
        # representation = "name"

        if representation == "electronic_structure" or representation == "periodic_table":
            arg = {"is_including_row": True}
            # arg = {}
            annot = False
        elif representation == "name":
            names = OFMHelper_structure2names(struct)
            arg = {"names": names}
            print(arg)
            annot = True

        ofmconverter = OFMConverter(representation=representation, **arg)
        ofm = ofmconverter.transform(struct)

        print("OFM is made.")

        OFMHelper_fig(ofm, action=["mean", "locals"], annot=annot)

    main()
