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

from .descriptor.RDF import RDFConverter

"""
requires 
spglib                             >=1.16.2
pymatgen                           >=2022.0.14
"""

sys.path.append("/home/kino/tmp/AkaiKKRPythonUtil/util")
if "test_script" not in sys.modules:
    from pyakaikkr import AkaikkrJob
    from pyakaikkr.Error import *
    from pyakaikkr import StructureSpeciesConverter


def _make_displc(anclr, displc=[0, 0, 0]):
    displc_list = []
    for anclr1 in anclr:
        displc1_list = []
        for z in anclr1:
            displc1_list.append(displc)
        displc_list.append(displc1_list)
    return displc_list


def _run_it(specx, ciffile, directory, structurefile="structure.json",
            displc=False, use_bravais=True, fmt="cif",
            inputcard="inputcard_geom", outputcard="out_geom.log"):
    """make inputcard from ciffile and run specx

    Args:
        specx (str): specx path
        ciffile (str): cif filename
        directory (str): directory to write inputcard and execute specx
        structurefile (str): json structure file. Defaults to "structure.json"
        displc (bool, optional): add displc in param or not. Defaults to False.
        use_bravis (bool, optional): use bravais lattice, Defaults to True. 
        fmt (bool, optional): ciffile format. Defaults to "cif".
        inputcard (str, optional): inputcard filename. Defaults to "inputcard_go".
        outputcard (str, optional): ouptutcard filename. Defaults to "out_go.log".

    Raises:
        CIF2KKRNoStructureError: failed to read structure
        KKRFailedExecutionError: failed to run specx, self.return_code is also made.

    Returns:
        bool: success or failed
    """
    os.makedirs(directory, exist_ok=True)

    meta = {"ciffile": ciffile, "directory": directory}
    filepath = os.path.join(directory, "meta.json")
    with open(filepath, "w") as f:
        json.dump(meta, f)

    job = AkaikkrJob(directory)
    param = job.default
    param["maxitr"] = 0  # stop after showing structures

    struc_param = None
    # first read as primitive. if it fails, read as conventional
    try:
        struc_param = job.read_structure(
            ciffile, use_bravais=use_bravais, cif_primitive=True, fmt=fmt)
    except CIF2KKRSpgDifferentError as err:
        print(err)
        print("WARNING: read structure again by cif_primitive=False")
        struc_param = job.read_structure(
            ciffile, use_bravais=use_bravais, cif_primitive=False, fmt=fmt)
    except CIF2KKRNsiteInconsistentError as err:
        print(err)
        print("WARNING: primitive and kkr nsite are inconsistent")
        try:
            struc_param = job.read_structure(
                ciffile, use_bravais=use_bravais, cif_primitive=False, fmt=fmt)
        except CIF2KKRNsiteInconsistentError as err:
            print(err)
            print("WARNING: primitive and kkr nsite are inconsistent")
            # possibility CIF2KKRNsiteInconsistentError occurs

    if struc_param is None:
        raise CIF2KKRNoStructureError("no structure is read")

    structurepath = os.path.join(directory, structurefile)
    with open(structurepath, "w") as f:
        json.dump(struc_param, f)
    print("save to", structurepath)

    param.update(struc_param)
    if displc:
        displc_param = {"displc": _make_displc(
            param["anclr"], displc=[0, 0, 0])}
        param.update(displc_param)

    param.update({"go": "geom"})
    job.make_inputcard(param, inputcard)
    print("saved to", os.path.join(directory, inputcard))

    if True:
        job.run(specx, inputcard, outputcard)
    else:
        cmd = f"cd {directory}; {specx} < {inputcard} > {outputcard}"
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            raise KKRFailedExecutionError("return_code={}".format(ret))

        stoppped_by_errtrp = job.check_stopped_by_errtrp(outputcard)
        if stoppped_by_errtrp:
            raise KKRFailedExecutionError("stopped by errtrp")

    return struc_param


def _read_output(directory: str, outputcard: str) -> Structure:
    """read outputcar in directory

    Args:
        directory (str): directory of job
        outputcard (str): outputcard filename

    Returns:
        pymatgen.core.Structure: structure
    """
    akaikkr = AkaikkrJob(directory)
    struc = akaikkr.make_pymatgenstructure(outputcard)
    return struc


if False:
    def make_rdf(structure: Structure, rcut=4.0) -> np.ndarray:
        """make radial distribution function

        Args:
            structure (pymatgen.core.structure): struture
            rcut (float, optional): cut off distance. Defaults to 4.0.

        Returns:
            np.ndarray: radial distribution function
        """
        nndistance = []
        for site in structure.sites:
            neighbor = structure.get_neighbors(site, r=rcut)
            for nn in neighbor:
                # nndistance.append(round(nn[1], nround))
                nndistance.append(nn[1])
        nndistance = np.array(nndistance)
        nndistance.sort()
        return nndistance


def show_equiv_matrix(prim_struc: Structure, input_analyzer: StructureMatcher) -> np.ndarray:
    """show qeuivalenet matrix

    Args:
        prim_struc (pymatgen.core.Structure): structure
        input_analyzer (pymatgen.analysis.structure_matcher.StructureMatcher): structure matcher

    Returns:
        np.ndarray: equivalenet matrix of sites by spacegroup operations
    """
    print("lattice", prim_struc.lattice)
    print("sites", prim_struc.sites)
    print("specis", prim_struc.species)

    species = prim_struc.species
    print("species", species)

    ops = input_analyzer.get_space_group_operations()
    print(ops)

    n = len(prim_struc.sites)
    # equiv_matrix = np.full((n, n), False)
    equiv_matrix = np.identity(n, dtype=bool)
    for i1 in range(n):
        for i2 in range(i1, n):
            site1 = PeriodicSite(
                species=species[i1], coords=prim_struc.sites[i1].frac_coords, lattice=prim_struc.lattice)
            site2 = PeriodicSite(
                species=species[i2], coords=prim_struc.sites[i2].frac_coords, lattice=prim_struc.lattice)

            eq = ops.are_symmetrically_equivalent([site1], [site2])
            equiv_matrix[i1, i2] = eq
    for i1 in range(n):
        for i2 in range(i1, n):
            equiv_matrix[i2, i1] = equiv_matrix[i1, i2]
    print(equiv_matrix)
    return equiv_matrix


def _make_both_rdf(input_rdf: np.ndarray, output_rdf: np.ndarray, bins: int) -> (np.ndarray, np.ndarray, np.ndarray, np.ndarray):
    """[summary]

    Args:
        input_rdf (np.ndarray): input RDF
        output_rdf (np.ndarray): output RDF
        bins (int): number of bins

    Returns:
        np.ndarray: input histogram
        np.ndarray: output histogram
        np.ndarray: input histogram bin edges
        np.ndarray: output histogram bin edges
    """
    if False:
        xmin = np.min([np.min(input_rdf), np.min(output_rdf)])
        xmax = np.max([np.max(input_rdf), np.max(output_rdf)])
        xrange = (xmin, xmax)
    else:
        xrange = []
        for op in [np.min, np.max]:
            xrange.append(op([op(input_rdf), op(output_rdf)]))
    input_hist, input_edges = np.histogram(
        input_rdf, bins=bins, range=xrange)
    output_hist, outputput_edges = np.histogram(
        input_rdf, bins=bins, range=xrange)
    return input_hist, output_hist, input_edges, outputput_edges


def _plot_rdf(input_rdf, output_rdf, directory, bins=200):
    """plot input and output RDF

    Args:
        input_rdf (nd.array): input RDF
        output_rdf (nd.array): output RDF
        directory (str): directory to save png
        bins (int, optional): number of bins. Defaults to 200.
    """
    if True:

        input_hist, output_hist, input_binedges, output_binedges = _make_both_rdf(
            input_rdf, output_rdf, bins)

        input_bin = (input_binedges[:-1] + input_binedges[1:])*0.5
        output_bin = (output_binedges[:-1] + output_binedges[1:])*0.5

        width = input_bin[1]-input_bin[0]
        fig, axes = plt.subplots(3, 1)
        ax = axes[0]
        ax.bar(input_bin, input_hist, width=width,
               color="blue")
        ax.set_title("cif (pymatgen)")
        ax = axes[1]
        ax.bar(output_bin, output_hist, width=width,
               color="blue")
        ax.set_title("kkr output")
        ax = axes[2]
        ax.axhline(0, linestyle="--")
        ax.bar(input_bin, output_hist-input_hist, width=width,)
        ax.set_xlabel("distance")
        ax.set_title("difference")

    else:
        fig, axes = plt.subplots(3, 1)

        ax = axes[0]
        ax.hist(input_rdf,
                color="blue", bins=bins)
        ax.set_title("cif (pymatgen)")

        ax = axes[1]
        ax.hist(output_rdf,
                color="red", bins=bins)
        ax.set_title("kkr output")

        ax = axes[2]
        ax.hist(input_rdf,
                alpha=0.5, color="blue", bins=bins)
        ax.hist(output_rdf,
                alpha=0.5, color="red", bins=bins)
        ax.set_xlabel("distance")
        ax.set_title("both")

    fig.tight_layout()
    filename = os.path.join(directory, "rdf.png")
    fig.savefig(filename)
    print("saved to", filename)
    fig.clf()
    plt.close(fig)


def _rdf_diff_is_zero(input_rdf, output_rdf, output_msg=False, bins=200):
    """whether the difference of RDF is zero or not

    Args:
        input_rdf (np.ndarray): input RDF
        output_rdf (np.ndarray): output RDF
        output_msg (str, optional): returns output message or bool. Defaults to False.
        bins (int, optional): number of bins. Defaults to 200.

    Returns:
        [type]: [description]
    """
    input_hist, output_hist, input_binedges, output_binedges = _make_both_rdf(
        input_rdf, output_rdf, bins)
    msg = []
    s = "rdf difference, bins= {}".format(bins)
    print(s)
    msg.append(s)
    diff = output_hist-input_hist
    print("_rdf_diff_is_zero", diff)
    msg.append(str(diff))

    if output_msg:
        return msg
    else:
        res = np.all(diff == 0)
        print("np.all(diff == 0)", res)
        return res


class CompareCifKkr:
    SUCCESS = "success"
    FAILED = "failed"
    NO_STRUCTURE = "no_structure"
    NO_ELEMENT = "unknown_element_in_input_file"
    FAILED_TO_GET_SPG_INFO = "failed_to_get_spginfo"
    FAILED_TO_RUN_SPECX = "failed_to_run_specx"
    NO_SPG_MATCH = "spg_not_match"

    def __init__(self, ciffile: str, specx: str, displc:  bool = False, fmt: str = "cif",
                 parent_directory: str = "output",
                 inputcard: str = "inputcard_geom", outputcard: str = "out_geom.log",
                 Vc: str = "Og", structurefile: str = "structure.json"):
        """Constructure

        Args:
            ciffile (str): a cif file
            specx (str): akaikkr program
            displc (bool, optional): include displc field or not. Defaults to False.
            fmt (str, optional): format of the cif file. Defaults to "cif".
            inputcard (str, optional): inputcard name. Defaults to "inputcard_geom".
            outputcard (str, optional): outputcard name. Defaults to "out_geom.log".
            Vc (str, optional): element treated as Z=0. Defaults to "Og".
            structurefile (str, optional): structure json file. Defaults to "structure.json".
        """
        self.ciffile = ciffile
        self.specx = specx
        self.displc = displc
        self.fmt = fmt
        self.Vc = Vc
        self.structurefile = structurefile

        self.inputcard = inputcard
        self.outputcard = outputcard

        self.all_done = False
        self.parent_directory = parent_directory
        hashstr = hashlib.md5(ciffile.encode()).hexdigest()
        self.directory = hashstr
        path = os.path.join(self.parent_directory, self.directory)
        os.makedirs(path, exist_ok=True)

        self.input_struc = None
        self.prim_stand_struc = None
        self.kkr_struc = None

        self.struc_param = None

    def matchedflag2msg(self, flag):
        """convert True|False to string message

        Args:
            flag (bool): structure is the same or not

        Returns:
            str: string message
        """
        if flag:
            return self.SUCCESS
        else:
            return self.FAILED

    def input_struc_to(self, fmt="cif", filename="input.cif"):
        filepath = os.path.join(self.parent_directory,
                                self.directory, filename)
        print("input_struc_to filepath", filepath)
        try:
            self.prim_stand_struc.to(fmt=fmt, filename=filepath)
        except IndexError:
            print(f"failed to write {filepath}\n"
                  "probably because of uknown element.\nbut continue.")

    def output_struc_to(self, fmt="cif", filename="output.cif"):
        filepath = os.path.join(self.parent_directory,
                                self.directory, filename)
        try:
            self.kkr_struc.to(fmt=fmt, filename=filepath)
        except IndexError:
            print(f"failed to write {filepath}\n"
                  "probably because of uknown element.\nbut continue.")

    def convert_and_compare(self, use_bravais=True, scale=False):
        """scale should be False

        Args:
            use_bravais (bool, optional): use bravais lattice. Defaults to True.
            scale (bool, optional): scale flag of StructureMatcher. Defaults to True.

        Raises:
            ValueError: [description]

        Returns:
            [type]: [description]
        """
        if self.fmt == "cif":
            try:
                parser = CifParser(self.ciffile)
            except AssertionError:
                self.msg = self.NO_STRUCTURE
                return self.NO_STRUCTURE

            try:
                self.conv_struc = parser.get_structures(primitive=False)[0]
            except ValueError:
                self.msg = self.NO_STRUCTURE
                return self.NO_STRUCTURE
            self.input_struc = self.conv_struc
            try:
                self.prim_struc = parser.get_structures(primitive=True)[0]
            except ValueError:
                self.msg = self.NO_STRUCTURE
                return self.NO_STRUCTURE
            self.input_struc = self.prim_struc

        elif self.fmt == "vasp" or self.fmt == "poscar":
            try:
                self.prim_struc = Structure.from_file(self.ciffile)
            except ValueError:
                self.msg = self.NO_STRUCTURE
                return self.NO_STRUCTURE
            self.conv_struc = self.prim_struc  # dummy
            self.input_struc = self.prim_struc

        self.conv_struc_spg_info = None
        try:
            self.conv_struc_spg_info = self.conv_struc.get_space_group_info()
        except TypeError:
            self.msg = "failed to get spginfo of conventional structure"
            return self.FAILED_TO_GET_SPG_INFO
        if self.conv_struc_spg_info is not None:
            print("get_structures(primitive=False), symmetry, nsite",
                  self.conv_struc_spg_info, len(self.conv_struc.sites))

        self.prim_struc_spg_info = None
        try:
            self.prim_struc_spg_info = self.prim_struc.get_space_group_info()
        except TypeError:
            self.msg = "failed to get spginfo of primitive structure"
            return self.FAILED_TO_GET_SPG_INFO
        if self.prim_struc_spg_info is not None:
            print("get_structures(primitive=True), symmetry, nsite",
                  self.prim_struc_spg_info, len(self.prim_struc.sites))

        analyzer = SpacegroupAnalyzer(self.conv_struc)
        try:
            self.prim_stand_struc = analyzer.get_primitive_standard_structure()
        except AttributeError as err:
            self.msg = "failed to use SpacegroupAnalyzer. Probably it is an alloy."
            print(self.msg)
            self.prim_stand_struc = None

        if self.prim_stand_struc is None:
            # try StructureSpeciesConverter
            structurespeciesconverter = StructureSpeciesConverter(
                self.conv_struc)
            substitutedstructure = structurespeciesconverter.structure
            analyzer = SpacegroupAnalyzer(
                substitutedstructure)
            cs_structure = analyzer.get_primitive_standard_structure()
            converted_structure = structurespeciesconverter.inverse_conversion(
                cs_structure)
            self.prim_stand_struc = converted_structure
            # delete local variables
            del substitutedstructure
            del cs_structure
            del analyzer

        self.prim_stand_struc_spg_info = self.prim_stand_struc.get_space_group_info()
        print("prim_stand_struc, symmetry, nsite",
              self.prim_stand_struc.get_space_group_info(),
              len(self.prim_stand_struc.sites))

        # use prim_stand_struc
        self.input_struc = self.prim_stand_struc
        self.input_struc_spg_info = self.prim_stand_struc_spg_info

        # @property self.input_struc = self.prim_stand_struc
        if self.input_struc is not None:
            self.input_struc_to()

        if self.input_struc:
            inputcard = self.inputcard
            outputcard = self.outputcard
            try:
                self.struc_param = _run_it(self.specx,  self.ciffile,
                                           directory=os.path.join(
                                               self.parent_directory, self.directory),
                                           structurefile=self.structurefile,
                                           displc=self.displc, use_bravais=use_bravais,
                                           fmt=self.fmt,
                                           inputcard=inputcard, outputcard=outputcard)
            except CIF2KKRGetStructureError as err:
                self.msg = str(err)
                return self.NO_STRUCTURE
            except CIF2KKRGetConventionalStandardStructureError as err:
                self.msg = str(err)
                return self.FAILED
            except CIF2KKRSpgDifferentError as err:
                self.msg = str(err)
                return self.FAILED
            except CIF2KKRCellShapeError as err:
                self.msg = str(err)
                return self.FAILED
            except CIF2KKRUnknownElementError as err:
                self.msg = str(err)
                return self.NO_ELEMENT
            except CIF2KKRNoStructureError as err:
                self.msg = str(err)
                return self.NO_STRUCTURE
            except KKRFailedExecutionError as err:
                self.msg = str(err)
                return self.FAILED_TO_RUN_SPECX

            try:
                self.kkr_struc = _read_output(os.path.join(
                    self.parent_directory, self.directory), outputcard)

            except ValueError:
                self.msg = self.NO_ELEMENT
                return self.NO_ELEMENT

            if self.kkr_struc:
                self.output_struc_to()
            else:
                self.msg = self.NO_STRUCTURE
                return self.NO_STRUCTURE

            self.kkr_struc_spg_info = self.kkr_struc.get_space_group_info()
            print("kkr_struc, symmetry, nsite",  self.kkr_struc_spg_info,
                  len(self.kkr_struc.sites))

        if self.input_struc and self.kkr_struc:
            matcher = StructureMatcher(scale=scale)
            mathced = False

            try:
                self.kkr_struc_spg_info = self.kkr_struc_spg_info
            except TypeError:
                print(
                    "failed to get space group symbol of structure. but continue.")
                raise TypeError

            print("kkr space group symbol, nsite",
                  self.kkr_struc_spg_info[1], len(self.kkr_struc.sites))

            input_struc = self.input_struc
            kkr_struc = self.kkr_struc

            print("input prim structure num_site",
                  input_struc.num_sites)
            print("kkr   prim structure num_site", kkr_struc.num_sites)

            # compare by size
            if len(input_struc.sites) != len(kkr_struc.sites):
                self.msg = "# of sites are different. skip showing diff of min_dist_matrix"
                matched = False
                return self.matchedflag2msg(matched)

            matched = matcher.fit(input_struc, kkr_struc)
            if matched:
                self.msg = "matched by structure matcher"
            else:
                # compare by rdf
                rdfconverter = RDFConverter()
                input_rdf = rdfconverter.make_nndistance(input_struc)['all']
                self.input_rdf = input_rdf
                output_rdf = rdfconverter.make_nndistance(kkr_struc)['all']
                self.output_rdf = output_rdf

                _plot_rdf(input_rdf, output_rdf, directory=os.path.join(
                    self.parent_directory, self.directory))
                if _rdf_diff_is_zero(input_rdf, output_rdf):
                    self.msg = "matched by rdf"
                    matched = True

            # optional check
            # spg may not match becaseu of a bug of spglib
            spgmatch = [True, True, True]
            if self.input_struc_spg_info[1] != self.kkr_struc_spg_info[1]:
                self.msg = "input and kkr space group is different {} != {}".format(
                    self.input_struc_spg_info[1], self.kkr_struc_spg_info[1])
                print("WARNING", self.msg)
                spgmatch[0] = False
            if self.prim_struc_spg_info[1] != self.kkr_struc_spg_info[1]:
                self.msg = "prim and kkr space group is different {} != {}".format(
                    self.prim_struc_spg_info[1], self.kkr_struc_spg_info[1])
                print("WARNING", self.msg)
                spgmatch[1] = False
            print(self.conv_struc_spg_info, self.kkr_struc_spg_info)
            if self.conv_struc_spg_info[1] != self.kkr_struc_spg_info[1]:
                self.msg = "conv and kkr space group is different {} != {}".format(
                    self.conv_struc_spg_info[1], self.kkr_struc_spg_info[1])
                print("WARNING", self.msg)
                spgmatch[2] = False
            spgmatch = np.array(spgmatch)

            if np.any(spgmatch == False):
                self.msg = self.matchedflag2msg(
                    matched)+", but "+self.NO_SPG_MATCH
            return self.matchedflag2msg(matched)

        self.msg = "structure isn't the same"
        return self.matchedflag2msg(matched)

    def log_msg(self, msg, filename=None):
        lines = []
        lines.append("ciffile = " + self.ciffile)
        lines.append("directory = " + self.directory)
        lines.append("msg = " + msg)

        try:
            lines.append("prim space group, number {}".format(
                self.prim_struc_spg_info))
        except AttributeError:
            pass
        try:
            lines.append("conv space group, number {}".format(
                self.conv_struc_spg_info))
        except AttributeError:
            pass
        try:
            lines.append("input  space group, number {}".format(
                self.input_struc_spg_info))
        except AttributeError:
            pass
        try:
            lines.append("output space group, number {}".format(
                self.kkr_struc_spg_info))
        except AttributeError:
            pass
        try:
            lines.append("input  structure len, {}".format(
                self.input_struc.num_sites))
        except AttributeError:
            pass
        try:
            lines.append("output structure len, {}".format(
                self.kkr_struc.num_sites))
        except AttributeError:
            lines.append("probably failed to make kkr structure")
        try:
            if self.input_struc.num_sites == self.kkr_struc.num_sites:
                lines.append("input  structure {}".format(self.input_struc))
                lines.append(("output structure {}".format(self.kkr_struc)))
        except AttributeError:
            pass

        if self.input_struc is not None and self.kkr_struc is not None:

            try:
                x = self.input_rdf
            except:
                rdfconverter = RDFConverter()
                self.input_rdf = rdfconverter.make_nndistance(self.input_struc)[
                    'all']

            try:
                x = self.output_rdf
            except:
                rdfconverter = RDFConverter()
                self.output_rdf = rdfconverter.make_nndistance(self.kkr_struc)[
                    'all']

            msg = _rdf_diff_is_zero(
                self.input_rdf, self.output_rdf, output_msg=True)
            lines.extend(msg)
            _plot_rdf(self.input_rdf, self.output_rdf, directory=os.path.join(
                self.parent_directory, self.directory))
        else:
            if self.input_struc is None:
                lines.append("input_struc is None")
            if self.kkr_struc is None:
                lines.append("kkr_struc is None")

        lines.append("")
        lines.append("")

        if filename is not None:
            with open(filename, "a") as f:
                f.write("\n".join(lines))
            print()
            print("appended to", filename)
            print()
        return lines

    def get_structure_param(self, remove_temperaryfiles=True):
        struc_file = os.path.join(self.parent_directory,
                                  self.directory, self.structurefile)
        with open(struc_file) as f:
            struc_param = json.load(f)

        if remove_temperaryfiles:
            path_remove = os.path.join(self.parent_directory,
                                       self.directory)
            shutil.rmtree(path_remove, ignore_errors=True)
            try:
                os.rmdir(self.parent_directory)
            except OSError as err:
                pass
        return struc_param
