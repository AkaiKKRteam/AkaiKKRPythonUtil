# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

#!/bin/env python

from .Error import *
from .Unit import *
import sys

import numpy as np
from pymatgen.io.cif import CifParser
from pymatgen.core import Structure, PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import Element
import pandas as pd

from .ElementKkr import ElementKKR

_kkr_bohr = Unit().length_au2ang


if False:
    try:
        from pymatgen.symmetry import kpath
        _use_kpath = False
    except ImportError:
        print("Warning: no kpath in pymatgen. kpath will be omitted.")
        _use_kpath = False


class _BravaisKKR:
    dict = {
        1: "trc",  # triclinic
        2: "trc",  # triclinic
        3: "sm",  # simple monoclinic
        4: "sm",  # simple monoclinic
        5: "bsm",  # base centered monoclinic
        6: "sm",  # simple monoclinic
        7: "sm",  # simple monoclinic
        8: "bsm",  # base centered monoclinic
        9: "bsm",  # base centered monoclinic
        10: "sm",  # simple monoclinic
        11: "sm",  # simple monoclinic
        12: "bsm",  # base centered monoclinic
        13: "sm",  # simple monoclinic
        14: "sm",  # simple monoclinic
        15: "bsm",  # base centered monoclinic
        16: "so",  # simple orthorhombic
        17: "so",  # simple orthorhombic
        18: "so",  # simple orthorhombic
        19: "so",  # simple orthorhombic
        20: "bso",  # base centered orthorhombic
        21: "bso",  # base centered orthorhombic
        22: "fco",  # face centered orthorhombic
        23: "bco",  # body centered orthorhombic
        24: "bco",  # body centered orthorhombic
        25: "so",  # simple orthorhombic
        26: "so",  # simple orthorhombic
        27: "so",  # simple orthorhombic
        28: "so",  # simple orthorhombic
        29: "so",  # simple orthorhombic
        30: "so",  # simple orthorhombic
        31: "so",  # simple orthorhombic
        32: "so",  # simple orthorhombic
        33: "so",  # simple orthorhombic
        34: "so",  # simple orthorhombic
        35: "bso",  # base centered orthorhombic
        36: "bso",  # base centered orthorhombic
        37: "bso",  # base centered orthorhombic
        38: "bso",  # base centered orthorhombic
        39: "bso",  # base centered orthorhombic
        40: "bso",  # base centered orthorhombic
        41: "bso",  # base centered orthorhombic
        42: "fco",  # face centered orthorhombic
        43: "fco",  # face centered orthorhombic
        44: "bco",  # body centered orthorhombic
        45: "bco",  # body centered orthorhombic
        46: "bco",  # body centered orthorhombic
        47: "so",  # simple orthorhombic
        48: "so",  # simple orthorhombic
        49: "so",  # simple orthorhombic
        50: "so",  # simple orthorhombic
        51: "so",  # simple orthorhombic
        52: "so",  # simple orthorhombic
        53: "so",  # simple orthorhombic
        54: "so",  # simple orthorhombic
        55: "so",  # simple orthorhombic
        56: "so",  # simple orthorhombic
        57: "so",  # simple orthorhombic
        58: "so",  # simple orthorhombic
        59: "so",  # simple orthorhombic
        60: "so",  # simple orthorhombic
        61: "so",  # simple orthorhombic
        62: "so",  # simple orthorhombic
        63: "bso",  # base centered orthorhombic
        64: "bso",  # base centered orthorhombic
        65: "bso",  # base centered orthorhombic
        66: "bso",  # base centered orthorhombic
        67: "bso",  # base centered orthorhombic
        68: "bso",  # base centered orthorhombic
        69: "fco",  # face centered orthorhombic
        70: "fco",  # face centered orthorhombic
        71: "bco",  # body centered orthorhombic
        72: "bco",  # body centered orthorhombic
        73: "bco",  # body centered orthorhombic
        74: "bco",  # body centered orthorhombic
        75: "st",  # simple tetragonal
        76: "st",  # simple tetragonal
        77: "st",  # simple tetragonal
        78: "st",  # simple tetragonal
        79: "bct",  # body centered tetragonal
        80: "bct",  # body centered tetragonal
        81: "st",  # simple tetragonal
        82: "bct",  # body centered tetragonal
        83: "st",  # simple tetragonal
        84: "st",  # simple tetragonal
        85: "st",  # simple tetragonal
        86: "st",  # simple tetragonal
        87: "bct",  # body centered tetragonal
        88: "bct",  # body centered tetragonal
        89: "st",  # simple tetragonal
        90: "st",  # simple tetragonal
        91: "st",  # simple tetragonal
        92: "st",  # simple tetragonal
        93: "st",  # simple tetragonal
        94: "st",  # simple tetragonal
        95: "st",  # simple tetragonal
        96: "st",  # simple tetragonal
        97: "bct",  # body centered tetragonal
        98: "bct",  # body centered tetragonal
        99: "st",  # simple tetragonal
        100: "st",  # simple tetragonal
        101: "st",  # simple tetragonal
        102: "st",  # simple tetragonal
        103: "st",  # simple tetragonal
        104: "st",  # simple tetragonal
        105: "st",  # simple tetragonal
        106: "st",  # simple tetragonal
        107: "bct",  # body centered tetragonal
        108: "bct",  # body centered tetragonal
        109: "bct",  # body centered tetragonal
        110: "bct",  # body centered tetragonal
        111: "st",  # simple tetragonal
        112: "st",  # simple tetragonal
        113: "st",  # simple tetragonal
        114: "st",  # simple tetragonal
        115: "st",  # simple tetragonal
        116: "st",  # simple tetragonal
        117: "st",  # simple tetragonal
        118: "st",  # simple tetragonal
        119: "bct",  # body centered tetragonal
        120: "bct",  # body centered tetragonal
        121: "bct",  # body centered tetragonal
        122: "bct",  # body centered tetragonal
        123: "st",  # simple tetragonal
        124: "st",  # simple tetragonal
        125: "st",  # simple tetragonal
        126: "st",  # simple tetragonal
        127: "st",  # simple tetragonal
        128: "st",  # simple tetragonal
        129: "st",  # simple tetragonal
        130: "st",  # simple tetragonal
        131: "st",  # simple tetragonal
        132: "st",  # simple tetragonal
        133: "st",  # simple tetragonal
        134: "st",  # simple tetragonal
        135: "st",  # simple tetragonal
        136: "st",  # simple tetragonal
        137: "st",  # simple tetragonal
        138: "st",  # simple tetragonal
        139: "bct",  # body centered tetragonal
        140: "bct",  # body centered tetragonal
        141: "bct",  # body centered tetragonal
        142: "bct",  # body centered tetragonal
        143: "hcp",  # hexagonal close packed
        144: "hcp",  # hexagonal close packed
        145: "hcp",  # hexagonal close packed
        146: "rhb",  # rhombohedral(trigonal)
        147: "hcp",  # hexagonal close packed
        148: "rhb",  # rhombohedral(trigonal)
        149: "hcp",  # hexagonal close packed
        150: "hcp",  # hexagonal close packed
        151: "hcp",  # hexagonal close packed
        152: "hcp",  # hexagonal close packed
        153: "hcp",  # hexagonal close packed
        154: "hcp",  # hexagonal close packed
        155: "rhb",  # rhombohedral(trigonal)
        156: "hcp",  # hexagonal close packed
        157: "hcp",  # hexagonal close packed
        158: "hcp",  # hexagonal close packed
        159: "hcp",  # hexagonal close packed
        160: "rhb",  # rhombohedral(trigonal)
        161: "rhb",  # rhombohedral(trigonal)
        162: "hcp",  # hexagonal close packed
        163: "hcp",  # hexagonal close packed
        164: "hcp",  # hexagonal close packed
        165: "hcp",  # hexagonal close packed
        166: "rhb",  # rhombohedral(trigonal)
        167: "rhb",  # rhombohedral(trigonal)
        168: "hcp",  # hexagonal close packed
        169: "hcp",  # hexagonal close packed
        170: "hcp",  # hexagonal close packed
        171: "hcp",  # hexagonal close packed
        172: "hcp",  # hexagonal close packed
        173: "hcp",  # hexagonal close packed
        174: "hcp",  # hexagonal close packed
        175: "hcp",  # hexagonal close packed
        176: "hcp",  # hexagonal close packed
        177: "hcp",  # hexagonal close packed
        178: "hcp",  # hexagonal close packed
        179: "hcp",  # hexagonal close packed
        180: "hcp",  # hexagonal close packed
        181: "hcp",  # hexagonal close packed
        182: "hcp",  # hexagonal close packed
        183: "hcp",  # hexagonal close packed
        184: "hcp",  # hexagonal close packed
        185: "hcp",  # hexagonal close packed
        186: "hcp",  # hexagonal close packed
        187: "hcp",  # hexagonal close packed
        188: "hcp",  # hexagonal close packed
        189: "hcp",  # hexagonal close packed
        190: "hcp",  # hexagonal close packed
        191: "hcp",  # hexagonal close packed
        192: "hcp",  # hexagonal close packed
        193: "hcp",  # hexagonal close packed
        194: "hcp",  # hexagonal close packed
        195: "sc",  # simple cubic
        196: "fcc",  # face centered cubic
        197: "bcc",  # body centered cubic
        198: "sc",  # simple cubic
        199: "bcc",  # body centered cubic
        200: "sc",  # simple cubic
        201: "sc",  # simple cubic
        202: "fcc",  # face centered cubic
        203: "fcc",  # face centered cubic
        204: "bcc",  # body centered cubic
        205: "sc",  # simple cubic
        206: "bcc",  # body centered cubic
        207: "sc",  # simple cubic
        208: "sc",  # simple cubic
        209: "fcc",  # face centered cubic
        210: "fcc",  # face centered cubic
        211: "bcc",  # body centered cubic
        212: "sc",  # simple cubic
        213: "sc",  # simple cubic
        214: "bcc",  # body centered cubic
        215: "sc",  # simple cubic
        216: "fcc",  # face centered cubic
        217: "bcc",  # body centered cubic
        218: "sc",  # simple cubic
        219: "fcc",  # face centered cubic
        220: "bcc",  # body centered cubic
        221: "sc",  # simple cubic
        222: "sc",  # simple cubic
        223: "sc",  # simple cubic
        224: "sc",  # simple cubic
        225: "fcc",  # face centered cubic
        226: "fcc",  # face centered cubic
        227: "fcc",  # face centered cubic
        228: "fcc",  # face centered cubic
        229: "bcc",  # body centered cubic
        230: "bcc",  # body centered cubic
    }

    @staticmethod
    def getType(group):
        if group in _BravaisKKR.dict:
            return _BravaisKKR.dict[group]
        else:
            return "aux"


class _TranslationKKR:
    matrix = {
        "sc":  np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
        "fcc": np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]),
        "bcc": np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
        "hcp": np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
        "rhb": np.array([[2.0, 1.0, 1.0], [-1.0, 1.0, 1.0], [-1.0, -2.0, 1.0]])/3.0,
        "st":  np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
        "bct": np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
        "so":  np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
        "fco": np.array([[0.0, 0.5, 0.5], [0.5, 0.0, 0.5], [0.5, 0.5, 0.0]]),
        "bco": np.array([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]]),
        "bso": np.array([[0.5, -0.5, 0.0], [0.5, +0.5, 0.0], [0.0, 0.0, +1.0]]),
        "sm":  np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
        "bsm": np.array([[0.5, -0.5, 0.0], [0.5, +0.5, 0.0], [0.0, 0.0, +1.0]]),
        "trc": np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
        "aux": np.array([[+1.0, 0.0, 0.0], [0.0, +1.0, 0.0], [0.0, 0.0, +1.0]]),
    }

    @staticmethod
    def getMatrix(brvtyp):
        if brvtyp in _TranslationKKR.matrix:
            return _TranslationKKR.matrix[brvtyp]
        else:
            return _TranslationKKR.matrix["aux"]


class _SiteKKR:
    def __init__(self, site, wyckoff, coords, Vc: str = "Og"):
        """initialization

        Args:
            site ([type]): [description]
            wyckoff ([type]): [description]
            coords ([type]): [description]
            Vc (str, optional): element treast as Z=0. Defaults to "Og".
        """
        def fold(c):
            c[0] = c[0] % 1.0
            c[1] = c[1] % 1.0
            c[2] = c[2] % 1.0
            if c[0] > 1.0-1.0e-8:
                c[0] = 0.0
            if c[1] > 1.0-1.0e-8:
                c[1] = 0.0
            if c[2] > 1.0-1.0e-8:
                c[2] = 0.0
            return c

        self.species = site.species
        self.frac_coords = fold(coords)

        conc_sum = 0.0  # sum of concentration.
        for element in site.species.elements:
            conc_sum += site.species.get_el_amt_dict()[element.symbol]

        if conc_sum < 1.0:
            self.name = site.species.reduced_formula + Vc + "_%s" % wyckoff
        else:
            self.name = site.species.reduced_formula + "_%s" % wyckoff

    def findSameCoords(self, vsite):
        def match_coords(ca, cb):
            return abs(ca[0]-cb[0]) < 1.0e-8 \
                and abs(ca[1]-cb[1]) < 1.0e-8 \
                and abs(ca[2]-cb[2]) < 1.0e-8

        found = False
        for site in vsite:
            if match_coords(self.frac_coords, site.frac_coords):
                found = True
                break

        return found

    def findSameName(self, vsite):
        found = False
        for site in vsite:
            if self.name == site.name:
                found = True
                break

        return found


def _get_spginfo_from_ciffile(filename):
    """extract space group data from a cif file
    returns None the cif file lacks it.

    Note:
        _symmetry_space_group_name_H-M, _symmetry_Int_Tables_number, _symmetry_cell_setting are obsolte, but are used widely. Therefore, they can be read.

    Args:
        filename (str): cif filename

    Returns:
        str: _space_group_name_H-M_alt
        int: _space_group_IT_number
        str: _space_group_crystal_system
    """
    with open(filename) as f:
        data = f.read().splitlines()
    name = None
    number = None
    cellsetting = None
    for line in data:
        s = line.split()
        if len(s) > 1:
            if s[0] == "_symmetry_space_group_name_H-M" or s[0] == "_space_group_name_H-M_alt":
                name = " ".join(s[1:])
            if s[0] == "_symmetry_Int_Tables_number" or s[0] == "_space_group_IT_number":
                number = int(s[1])
            if s[0] == "_symmetry_cell_setting" or s[0] == "_space_group_crystal_system":
                cellsetting = s[1]
    # debug print
    if False:
        print("ciffile: name, number, cellsetting", name, number, cellsetting)
    return name, number, cellsetting


class StructureSpeciesConverter:
    def __init__(self, structure):
        # make unit species dictionary
        species_dict = []
        for site in structure.sites:
            newitem = dict(site.species.as_dict())
            if newitem not in species_dict:
                species_dict.append(newitem)

        # make conversion table
        if True:
            conv_table = {}
            Z = 0
            for content in species_dict:
                Z += 1
                if Z == 2 or Z == 10 or Z == 18 or Z == 36 or Z == 54 or Z == 86:
                    # skip noble gas
                    Z += 1
                if Z > 103:
                    raise CIF2KKRTooManyTypesError("too many type definitions")
                elm = Element("H").from_Z(Z)
                elm_str = str(elm)
                conv_table[elm_str] = content
        else:
            conv_table = {}
            for Z, content in zip(range(1, 104), species_dict):
                elm = Element("H").from_Z(Z)
                elm_str = str(elm)
                conv_table[elm_str] = content

        self.conv_table = conv_table

        # make a dummay struture
        new_species = []
        for site in structure.sites:
            newitem = dict(site.species.as_dict())
            keys = [k for k, v in conv_table.items() if v == newitem]
            if len(keys) != 1:
                print("strange keys", keys)
                raise ValueError
            new_species.append(keys[0])

        # make fractional coordinates
        frac_coords = []
        for site in structure.sites:
            frac_coords.append(site.frac_coords)

        new_structure = Structure(
            lattice=structure.lattice, species=new_species, coords=frac_coords)
        self.substituted_structure = new_structure

    @ property
    def structure(self):
        return self.substituted_structure

    def inverse_conversion(self, structure):
        conv_table = self.conv_table
        # make a dummay struture
        new_species = []
        for site in structure.sites:
            newitem = dict(site.species.as_dict())
            if len(newitem.keys()) != 1:
                print("keys >1", newitem)
                raise ValueError

            key = list(newitem.keys())[0]

            if newitem[key] != 1.0:
                raise ValueError("value error in inverse_conversion")
            new_species.append(conv_table[key])

        # make fractional coordinates
        frac_coords = []
        for site in structure.sites:
            frac_coords.append(site.frac_coords)

        new_structure = Structure(
            lattice=structure.lattice, species=new_species, coords=frac_coords)
        return new_structure


def _show_equiv_matrix(structure, input_analyzer, wy):
    """obsolete"""
    species = structure.species
    ops = input_analyzer.get_space_group_operations()

    n = len(structure.sites)
    # equiv_matrix = np.full((n, n), False)
    equiv_matrix = np.identity(n, dtype=bool)
    for i1 in range(n):
        for i2 in range(i1, n):
            site1 = PeriodicSite(
                species=species[i1], coords=structure.sites[i1].frac_coords, lattice=structure.lattice)
            site2 = PeriodicSite(
                species=species[i2], coords=structure.sites[i2].frac_coords, lattice=structure.lattice)

            eq = ops.are_symmetrically_equivalent([site1], [site2])
            equiv_matrix[i1, i2] = eq
    for i1 in range(n):
        for i2 in range(i1, n):
            equiv_matrix[i2, i1] = equiv_matrix[i1, i2]
    # make indeces and columns
    namelist = []
    for specie, wy in zip(species, wy):
        namelist.append("{},{}".format(specie, wy))
    df = pd.DataFrame(equiv_matrix, index=namelist, columns=namelist)
    print(df)

    uniq_name_wy_list = list(set(namelist))
    print(uniq_name_wy_list)
    uniq_name_wy_count = {}
    for key in uniq_name_wy_list:
        uniq_name_wy_count[key] = 0

    # make submatrix step by step
    nlist = list(range(len(namelist)))
    for i in nlist:
        flaglist = df.iloc[i, :].values
        ilist = np.where(flaglist == True)[0]
        dfs = df.iloc[ilist, ilist]
        name_wy = dfs.columns[0]
        s = name_wy.split(",")
        if len(s) == 2:
            print(dfs)
            uniq_name_wy_count[name_wy] += 1
            id_ = uniq_name_wy_count[name_wy]
            new_name_wy = "{},{}".format(name_wy, id_)
            for j in ilist:
                namelist[j] = new_name_wy
            df.index = namelist
            df.columns = namelist

    print(df)
    elm_list = []
    wy_list = []
    id_list = []
    for name in df.columns:
        s = name.split(",")
        elm_list.append(s[0])
        wy_list.append(s[1])
        id_list.append(s[2])

    return elm_list, wy_list, id_list


def _get_uniq_wyckoff(structure):
    analyzer = SpacegroupAnalyzer(structure)

    wyckoffs = analyzer.get_symmetry_dataset()["wyckoffs"]
    equiv_atoms = analyzer.get_symmetry_dataset()["equivalent_atoms"]

    if len(wyckoffs) != len(equiv_atoms) or len(wyckoffs) != len(structure.sites):
        print("len(wyckoffs)={}, ".format(len(wyckoffs)) +
              "len(equiv_atoms)={}, ".format(len(equiv_atoms)) +
              "len(structure.sites)={}".format(len(structure.sites)))
        raise ValueError(
            "len(wyckoffs)!=len(equiv_atoms), possibly !=len(structure.sites)")

    wyckoffs_conv = []
    for site, wy, eq in zip(structure.sites, wyckoffs, equiv_atoms):
        mul = np.count_nonzero(equiv_atoms == eq)
        wyckoffs_conv.append("{}{}_{}".format(mul, wy, str(eq)))

    # debug print
    if False:
        print(wyckoffs_conv)
        for site, wy in zip(structure.sites, wyckoffs_conv):
            print("spceie, wy,eq", site.as_dict()["species"], wy)

    return wyckoffs_conv

    # obsolete algorithm
    # Kino keeps it for future use
    if True:
        speciesconverter = StructureSpeciesConverter(structure)
        substitutedstructure = speciesconverter.structure
    else:
        substitutedstructure = structure
    print(substitutedstructure)

    elm_list, wy_list, id_list = _show_equiv_matrix(
        substitutedstructure, analyzer, wyckoffs)
    print(elm_list, wy_list, id_list)

    converted_structure = speciesconverter.inverse_conversion(
        substitutedstructure)

    print(converted_structure)

    wyckoff_conv = []
    for wy, id_ in zip(wy_list, id_list):
        wyckoff_conv.append("{}_{}".format(wy, id_))

    for site, wy in zip(converted_structure.sites, wyckoff_conv):
        namedict = site.species.as_dict()
        name = None
        if len(namedict.keys()) == 1:
            key = list(namedict.keys())[0]
            if namedict[key] == 1.0:
                name = key
        if name is None:
            name = str(site.species)
        print(str(site.species), "_".join([name, wy]))

    return wyckoff_conv


def _found_unknown_elements(structure):
    elementkkr = ElementKKR(Vc=None)
    elementkkr = list(elementkkr.dict.keys())

    sites = structure.sites
    for site in sites:
        for elm in site.species.elements:
            elm = str(elm)
            if elm not in elementkkr:
                print("unknown element", elm)
                print("known elements are", elementkkr)
                return True
    return False


def _show_cell_parameters(structure_conv):
    print("# conventional cell")
    print("   a=%9.5f b=%9.5f c=%9.5f(a.u)" %
          (structure_conv.lattice.a,
           structure_conv.lattice.b,
           structure_conv.lattice.c))
    print("   alpha=%9.5f beta=%9.5f gamma=%9.5f(degree)" %
          (structure_conv.lattice.alpha,
           structure_conv.lattice.beta,
           structure_conv.lattice.gamma))


def _show_lattice_parameters(lattice_constant, structure_conv, lattice_prim):
    print("# lattice constant a %9.5f [angstrom]" % lattice_constant)
    print("# conventional translation vectors (in units of a)")
    print("   a=(%9.5f%9.5f%9.5f)" %
          (structure_conv.lattice.matrix[0][0]/lattice_constant,
           structure_conv.lattice.matrix[0][1]/lattice_constant,
           structure_conv.lattice.matrix[0][2]/lattice_constant))
    print("   b=(%9.5f%9.5f%9.5f)" %
          (structure_conv.lattice.matrix[1][0]/lattice_constant,
           structure_conv.lattice.matrix[1][1]/lattice_constant,
           structure_conv.lattice.matrix[1][2]/lattice_constant))
    print("   c=(%9.5f%9.5f%9.5f)" %
          (structure_conv.lattice.matrix[2][0]/lattice_constant,
           structure_conv.lattice.matrix[2][1]/lattice_constant,
           structure_conv.lattice.matrix[2][2]/lattice_constant))
    volume = np.dot(np.cross(structure_conv.lattice.matrix[0],
                             structure_conv.lattice.matrix[1]),
                    structure_conv.lattice.matrix[2])
    print("   volume=   %10.5f(a.u.)" % (volume/_kkr_bohr**3))

    print("# primitive translation vectors (in units of a)")
    print("   a=(%9.5f%9.5f%9.5f)" %
          (lattice_prim[0][0]/lattice_constant,
           lattice_prim[0][1]/lattice_constant,
           lattice_prim[0][2]/lattice_constant))
    print("   b=(%9.5f%9.5f%9.5f)" %
          (lattice_prim[1][0]/lattice_constant,
           lattice_prim[1][1]/lattice_constant,
           lattice_prim[1][2]/lattice_constant))
    print("   c=(%9.5f%9.5f%9.5f)" %
          (lattice_prim[2][0]/lattice_constant,
           lattice_prim[2][1]/lattice_constant,
           lattice_prim[2][2]/lattice_constant))

    volume = np.dot(
        np.cross(lattice_prim[0], lattice_prim[1]), lattice_prim[2])
    print("   volume=   %10.5f(a.u.)" % (volume/_kkr_bohr**3))


def _show_atomic_position(sites_prim, lattice_prim, lattice_constant):
    print("# atomic positions (in units of a)  %d atoms" % len(sites_prim))
    for site in sites_prim:
        position = np.dot(site.frac_coords, lattice_prim)/lattice_constant
        if abs(position[0]) < 1e-6:
            position[0] = 0.0
        if abs(position[1]) < 1e-6:
            position[1] = 0.0
        if abs(position[2]) < 1e-6:
            position[2] = 0.0
        print("   position=%13.8f%13.8f%13.8f  type=%s" %
              (position[0], position[1], position[2], site.name))


def ak_cif2kkrparam(filename: str, use_bravais: bool = True, use_primitive: bool = True,
                    cif_primitive: bool = True,
                    fmt: str = "cif", Vc: str = "Og",
                    show_detail: bool = False):
    """
        check whether the cif space group number is the same as that of spglib
        check whether the number of sites of the cif file is the same as that of this routine

        If use_bravais is True, use_primitive is set to True

        if fmt is "cif", CifParser is used. In the other fmt, Structure.from_file() is used.


    Args:
        filename (str): cif filename
        use_bravais (bool, optional): use bravias lattice. Defaults to True.
        use_primitive (bool, optional): use primitive cell. Defaults to True.
        cif_primitive (bool, optional): read the cif file as primitive cell. Defaults to True.
        fmt (str, optional): filename format. Defaults to "cif".
        show_detail (bool, optional): [description]. Defaults to False.

    Raises:
        CIF2KKRGetStructureError: failed to read structure via CifParser
        CIF2KKRUnknownElementError: unknown element in the cif file
    """

    if fmt == "cif":
        parser = CifParser(filename)
        try:
            print("cif_primitive=", cif_primitive)
            structure_work = parser.get_structures(primitive=cif_primitive)[0]
        except ValueError:
            raise CIF2KKRGetStructureError("failed in parser.get_structures.\n"
                                           + "please correct occupancies and coordinates.")
    else:
        try:
            structure_work = Structure.from_file(filename)
        except ValueError:
            raise CIF2KKRGetStructureError("failed in Struture.from_file.\n"
                                           + "please check occupancies and coordinates.")

    if _found_unknown_elements(structure_work):
        raise CIF2KKRUnknownElementError("unknown element in the cif file")

    analyzer = SpacegroupAnalyzer(structure_work)
    # analyzer = SpacegroupAnalyzer(structure_work,symprec=0.001, angle_tolerance=0.5) # bad result

    try:
        structure_conv = analyzer.get_conventional_standard_structure()
    except TypeError:
        raise CIF2KKRGetConventionalStandardStructureError("failed in analyzer.get_conventional_standard_structure.\n"
                                                           + "please correct occupancies and coordinates.")

    param = {}
    if use_bravais:
        use_primitive = True

    # setup primitive cell vectors.
    if use_primitive:

        if fmt == "cif":
            _, cif_number, _ = _get_spginfo_from_ciffile(filename)
        if False:
            spginfo = structure_work.get_space_group_info()
            number = spginfo[1]
        else:
            number = analyzer.get_space_group_number()
        if fmt == "cif":
            print("cif symmetry, spg lib symmetry", cif_number, number)
            if cif_number != number:
                print("WARNING: spg number in the cif file != spg number from spglib")
                # raise CIF2KKRSpgDifferentError(
                #    "spg number in the cif file != spg number from spglib")

        brvtyp = _BravaisKKR.getType(number)

        if show_detail:
            print("# space group")
            print("   number=%d  bravais=%s kkr_brvtyp=%s" % (
                number,
                analyzer.get_crystal_system(), brvtyp))

        matrix = _TranslationKKR.getMatrix(brvtyp)
        lattice_prim = np.dot(matrix, structure_conv.lattice.matrix)

        if use_bravais:
            param["brvtyp"] = brvtyp
        else:
            param["brvtyp"] = "aux"
    else:
        lattice_prim = structure_conv.lattice.matrix  # overwrite by conventional cell
        param["brvtyp"] = "aux"

    del structure_work
    del analyzer

    if show_detail:
        _show_cell_parameters(structure_conv)

    inv_lattice_prim = np.linalg.inv(lattice_prim)

    lattice_constant = structure_conv.lattice.a

    param["a"] = lattice_constant/_kkr_bohr
    param["c/a"] = structure_conv.lattice.c/lattice_constant
    param["b/a"] = structure_conv.lattice.b/lattice_constant
    param["alpha"] = structure_conv.lattice.alpha
    param["beta"] = structure_conv.lattice.beta
    param["gamma"] = structure_conv.lattice.gamma

    # change lattice constant A and angle alpha from hex to rhb.
    if param["brvtyp"] == "rhb":
        aa = np.dot(lattice_prim[0], lattice_prim[0])
        ab = np.dot(lattice_prim[0], lattice_prim[1])
        lattice_constant = np.sqrt(aa)
        theta = np.arccos(ab/aa)

        param["a"] = lattice_constant/_kkr_bohr
        param["c/a"] = 0.0  # not used in KKR
        param["b/a"] = 0.0  # not used in KKR
        param["alpha"] = theta*180.0/np.pi
        param["beta"] = 0.0  # not used in KKR
        param["gamma"] = 0.0  # not used in KKR

    # check lattice lengthes and angles.
    if param["brvtyp"] in ["sc", "fcc", "bcc"]:
        if structure_conv.lattice.a == structure_conv.lattice.b and \
           structure_conv.lattice.b == structure_conv.lattice.c and \
           abs(structure_conv.lattice.alpha-90.0) < 1e-8 and \
           abs(structure_conv.lattice.beta - 90.0) < 1e-8 and \
           abs(structure_conv.lattice.gamma-90.0) < 1e-8:
            pass
        else:
            raise CIF2KKRCellShapeError("invalid cell shape. "
                                        + "please check cell lengthes and angles.")
    elif param["brvtyp"] in ["hcp", "rhb"]:
        if structure_conv.lattice.a == structure_conv.lattice.b and \
           abs(structure_conv.lattice.alpha - 90.0) < 1e-8 and \
           abs(structure_conv.lattice.beta - 90.0) < 1e-8 and \
           abs(structure_conv.lattice.gamma-120.0) < 1e-8:
            pass
        else:
            raise CIF2KKRCellShapeError("invalid cell shape. "
                                        + "please check cell lengthes and angles.")
    elif param["brvtyp"] in ["st", "bct", "so", "fco", "bco", "bso"]:
        if abs(structure_conv.lattice.alpha-90.0) < 1e-8 and \
           abs(structure_conv.lattice.beta - 90.0) < 1e-8 and \
           abs(structure_conv.lattice.gamma-90.0) < 1e-8:
            pass
        else:
            raise CIF2KKRCellShapeError("invalid cell shape. "
                                        + "please check cell angles.")
    elif param["brvtyp"] in ["sm", "bsm"]:
        if abs(structure_conv.lattice.alpha-90.0) < 1e-8 and \
           abs(structure_conv.lattice.gamma-90.0) < 1e-8:
            pass
        else:
            raise CIF2KKRCellShapeError("invalid cell shape. "
                                        + "please check cell angles.")
    elif param["brvtyp"] in ["trc"]:
        pass

    # get rotation matrix which rorates
    # the convensional a vector in x-axis,
    # the conventional b vector on x-y plane.
    # rotate the primitive vectors by the rotation matrix.
    if param["brvtyp"] in ["sm", "bsm", "bso", "trc"]:
        va = structure_conv.lattice.matrix[0]  # A-axis
        vb = structure_conv.lattice.matrix[1]  # B-axis
        vc = structure_conv.lattice.matrix[2]  # C-axis

        ea = va/np.sqrt(np.dot(va, va))
        eb = vb/np.sqrt(np.dot(vb, vb))
        ec = np.cross(ea, eb)
        ec = ec/np.sqrt(np.dot(ec, ec))
        eb = np.cross(ec, ea)
        eb = eb/np.sqrt(np.dot(eb, eb))

        rot_mat = np.array([ea, eb, ec]).T

        lattice_prim = np.dot(lattice_prim, rot_mat)

    param["r1"] = lattice_prim[0]/lattice_constant
    param["r2"] = lattice_prim[1]/lattice_constant
    param["r3"] = lattice_prim[2]/lattice_constant
    # change np.ndarray to list to treat as json
    for key in ["r1", "r2", "r3"]:
        param[key] = param[key].tolist()

    if show_detail:
        _show_lattice_parameters(
            lattice_constant, structure_conv, lattice_prim)

    wyckoff_conv = _get_uniq_wyckoff(structure_conv)

    # construct primitive sites.
    sites_prim = []
    for site, wyckoff in zip(structure_conv.sites, wyckoff_conv):
        position = np.dot(site.frac_coords, structure_conv.lattice.matrix)
        coords = np.dot(position, inv_lattice_prim)
        site_prim = _SiteKKR(site, wyckoff, coords)

        if not site_prim.findSameCoords(sites_prim):
            sites_prim.append(site_prim)

    # choose unique site types.
    sites_type = []
    for site_prim in sites_prim:
        if not site_prim.findSameName(sites_type):
            sites_type.append(site_prim)

    # calculate properties of each site type.
    param["ntyp"] = len(sites_type)
    param["type"] = []
    param["ncmp"] = []
    param["rmt"] = []
    param["field"] = []
    param["mxl"] = []
    param["anclr"] = []
    param["conc"] = []

    elementkkr = ElementKKR(Vc=Vc)
    for site in sites_type:
        rmt_avg = 0.0  # averaged radius of maffin-tin potential [lattice a].
        field_max = 0.0  # maximum magnetic field[Ry].
        lmax_max = 0   # maximum angular momentum.
        conc_sum = 0.0  # sum of concentration.

        for element in site.species.elements:
            rmt_avg += elementkkr.getAtomicRMT(element.symbol)
            field = elementkkr.getAtomicField(element.symbol)
            lmax = elementkkr.getAtomicLMax(element.symbol)
            conc_sum += site.species.get_el_amt_dict()[element.symbol]
            if field_max < abs(field):
                field_max = abs(field)
            if lmax_max < lmax:
                lmax_max = lmax
        rmt_avg /= len(site.species.elements)

        param["type"].append(site.name)
        if conc_sum < 1.0:
            # add vacancy as element Vc
            param["ncmp"].append(len(site.species.elements)+1)
        else:
            param["ncmp"].append(len(site.species.elements))

        param["rmt"].append(rmt_avg)
        param["field"].append(field_max)
        param["mxl"].append(lmax_max)

        anclr = []
        conc = []

        for element in site.species.elements:
            anclr.append(elementkkr.getAtomicNumber(element.symbol))
            conc.append(site.species.get_el_amt_dict()[element.symbol]*100.0)

        if conc_sum < 1.0:
            # add vacancy as element Vc
            anclr.append(elementkkr.getAtomicNumber(Vc))
            conc.append((1.0-conc_sum)*100.0)

        param["anclr"].append(anclr)
        param["conc"].append(conc)

    # calculate atomic positions.
    param["natm"] = len(sites_prim)
    param["atmicx"] = []

    for site in sites_prim:
        param["atmicx"].append([
            "%10.8fa" % site.frac_coords[0],
            "%10.8fb" % site.frac_coords[1],
            "%10.8fc" % site.frac_coords[2],
            site.name])

    if show_detail:
        _show_atomic_position(sites_prim, lattice_prim, lattice_constant)

    nsite_kkr = len(param["atmicx"])

    def _get_nsite_prim(filename: str, fmt: str):
        """get nsite of input and primitive standard structure"""
        if fmt == "cif":
            print("conventional standard structure")
            parser = CifParser(filename)
            structure = parser.get_structures(primitive=True)[0]
            cif_nsite = len(structure.sites)

            structure = parser.get_structures(primitive=False)[0]
        else:
            structure = Structure.from_file(filename)
            cif_nsite = len(structure.sites)

        speciesconverter = StructureSpeciesConverter(structure)
        substitutedstructure = speciesconverter.structure
        analyzer = SpacegroupAnalyzer(substitutedstructure)
        prim_stand_structure = analyzer.get_primitive_standard_structure()
        prim_stand_nsite = len(prim_stand_structure)
        return cif_nsite, prim_stand_nsite

    cif_nsite, prim_stand_nsite = _get_nsite_prim(filename, fmt=fmt)
    print("prim_stand_nsite={}, nsite_kkr={}".format(prim_stand_nsite, nsite_kkr))
    if nsite_kkr != prim_stand_nsite:
        raise CIF2KKRNsiteInconsistentError(
            "prim_stand_nsite={} != nsite_kkr={}".format(prim_stand_nsite, nsite_kkr))

    # calculate bandmap k-path for primitive cell.
    # Kino keeps it for future use
    if False:
        species_prim = []
        coords_prim = []
        for site in sites_prim:
            species_prim.append(site.species)
            coords_prim.append(site.frac_coords)
        # end for site
        structure_prim = Structure(lattice_prim, species_prim, coords_prim)
        try:
            kpseek = kpath.KPathSeek(structure_prim)

            kpath_lines = []
            knpoints = 150
            kpath_lines.append("%d" % knpoints)

            kpath_line = "c ---"
            # only the first succesive segments
            for ksymbol in kpseek.kpath["path"][0]:
                kpath_line += " %s" % ksymbol
            # end for ksymbol
            kpath_lines.append(kpath_line)

            # only the first succesive segments
            for ksymbol in kpseek.kpath["path"][0]:
                kpath_lines.append("c --- %s" % ksymbol)
                kvector = kpseek.kpath["kpoints"][ksymbol]
                kpath_lines.append("%.2f %.2f %.2f" %
                                   (kvector[0], kvector[1], kvector[2]))
            # end for ksymbol

            param["kpath_raw"] = kpath_lines
        except TypeError:
            print("warning: failed in kpath.KPathSeek.")
#            raise ValueError( "failed in kpath.KPathSeek." )

    return param
