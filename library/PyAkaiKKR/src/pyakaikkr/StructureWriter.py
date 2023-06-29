# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


def _make_spcecies_shortname(label: str, species: list):
    """make secies short name

    Args:
        label (str): label
        species (list): a list of species

    Returns:
        str: shortname
    """
    if len(species) == 1:
        return label
    else:
        _label = []
        for specie in species:
            elm = specie["element"]
            # occu = str(int(specie["occu"]*100))
            occu = specie["occu"]
            _label.append("{}{:.02f}".format(elm, occu))
        return "".join(_label)


def _kkr2cif(structure, lattice="direct", fmt="cif",
             outputfilename="structure.cif", include_Vc=True):
    """make a cif file.

    if outputfilename == None, print to stdout

    Args:
        struc (Structure): pymatgen.core.Structure
        lattice (str, optional): lattice conversion. Defaults to "direct".
        outputfilename (str, optional): output filename. Defaults to "crystal.cif".
        include_Vc (bool, optional): inlucde Z=0 sites. Defaults to True

    """
    lines = []

    chemicalformula = structure.composition.formula
    msg = "data_{}".format(chemicalformula.replace(" ", ""))
    lines.append(msg)

    lines.append("_symmetry_space_group_name_H-M   'P 1'")

    for abckey, abcvalue in zip("abc", [structure.lattice.a, structure.lattice.b, structure.lattice.c]):
        lines.append("_cell_length_{} {}".format(abckey, abcvalue))

    for abckey, abcvalue in zip(["alpha", "beta", "gamma"],
                                [structure.lattice.alpha, structure.lattice.beta,
                                 structure.lattice.gamma]):
        lines.append("_cell_angle_{} {}".format(abckey, abcvalue))

    lines.append("_chemical_formula_structural {}".format(
        structure.composition.formula))
    lines.append("_chemical_formula_sum {}".format(
        structure.composition.reduced_formula))
    lines.append("_cell_formula_units_Z 1")

    msg = """loop_
 _symmetry_equiv_pos_site_id
 _symmetry_equiv_pos_as_xyz
  1  'x, y, z'

loop_
 _atom_site_type_symbol
 _atom_site_label
 _atom_site_symmetry_multiplicity
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 _atom_site_occupancy
"""

    for s in msg.split("\n"):
        lines.append(s)

    iserial = 0
    for site in structure.sites:
        frac_coords = site.frac_coords.tolist()
        frac_coords = list(map(str, frac_coords))
        species_dict = dict(site.species.as_dict())
        if len(species_dict.keys()) == 0 and include_Vc:
            elm = "Vc"
            occup = 1.0
            line = " {} {}{} 1 {} {}".format(
                elm, elm, iserial, " ".join(frac_coords), occup)
            iserial += 1
            lines.append(line)
        else:
            for elm, occup in site.species.as_dict().items():
                line = " {} {}{} 1 {} {}".format(
                    elm, elm, iserial, " ".join(frac_coords), occup)
                iserial += 1
                lines.append(line)
    lines.append("")

    if outputfilename is not None:
        with open(outputfilename, "w") as f:
            f.write("\n".join(lines))
        print("saved to", outputfilename)
    else:
        print("\n".join(lines))


def _writer(struc: Structure, lattice="direct", fmt="cif",
            outputfilename="structure.cif"):
    """make a cif|poscar file.
        lattice=primitive|conventional|direct

        lattice=primitive generates primitive standard structure.
        lattice=conventional generates conventional standard structure.
        lattice=direct outputs the structure in struc.

    Args:
        struc (Structure): pymatgen.core.Structure
        lattice (str, optional): lattice conversion. Defaults to "primtitive".
        outputfilename (str, optional): output filename. Defaults to "crystal.cif".

    Raises:
        ValueError: unknown lattice
    """
    if lattice == "direct":
        _struc = struc
    elif lattice == "primtitive":
        sga = SpacegroupAnalyzer(struc)
        _struc = sga.get_primitive_standard_structure()
    elif lattice == "convential":
        sga = SpacegroupAnalyzer(struc)
        _struc = sga.get_conventional_standard_structure()
    else:
        raise ValueError("unknown lattice={}".format(lattice))

    if outputfilename is not None:
        _struc.to(fmt=fmt, filename=outputfilename)
        print("saved to", outputfilename)
    else:
        print("console mode not supported in the poscar format.")


def _writer_xsf(struc: Structure, lattice="standard", outputfilename="structure.xsf"):
    """make a xsf file.
    lattice = standard|direct

    lattice='standard' generates primitive and conventional matrices.
    lattice='direct' make primitive and conventional matrices from the given struc.

    Args:
        struc (Structure): pymatgen.core.Structure
        lattice (str, optional): lattice type. Defaults to "standard".
        outputfilename (str, optional): output filename. Defaults to "crystal.xsf".

    Raises:
        ValueError: unknown lattice
    """
    if lattice == "direct":
        prim_struc_dic = struc.as_dict()
        conv_struc_dic = struc.as_dict()
    elif lattice == "standard":
        sga = SpacegroupAnalyzer(struc)
        prim_struc = sga.get_primitive_standard_structure()
        prim_struc_dic = prim_struc.as_dict()
        if False:
            for dkey, content in prim_struc_dic.items():
                if dkey == "lattice":
                    print(dkey, content["matrix"])
                if dkey == "sites":
                    for sites in content:
                        print(sites)
        conv_struc = sga.get_conventional_standard_structure()
        conv_struc_dic = conv_struc.as_dict()
    else:
        raise ValueError("unknown lattice={}".format((lattice)))

    lines = []
    lines = ["CRYSTAL"]
    lines.append("PRIMVEC")
    for lattice in prim_struc_dic["lattice"]["matrix"]:
        lattice = list(map(str, lattice))
        lines.append(" ".join(lattice))
    lines.append("CONVVEC")
    for lattice in conv_struc_dic["lattice"]["matrix"]:
        lattice = list(map(str, lattice))
        lines.append(" ".join(lattice))
    lines.append("PRIMCOORD")
    nsite = len(prim_struc_dic["sites"])
    lines.append(" ".join([str(nsite), "1"]))
    for sites in prim_struc_dic["sites"]:
        label = _make_spcecies_shortname(sites["label"], sites["species"])
        xyz = list(map(str, sites["xyz"]))
        s = [label]
        s.extend(xyz)
        lines.append(" ".join(s))

    if outputfilename is not None:
        with open(outputfilename, "w") as f:
            f.write("\n".join(lines))
            print("saved to", outputfilename)
    else:
        print("\n".join(lines))


class StructureWriter:
    def __init__(self, structure):
        self.structure = structure

    def cif(self, latticetype="direct", outputfilename="structure.cif"):
        """write as cif file

        Args:
            latticetype (str, optional): lattice type. Defaults to "direct".
            outputfilename (str, optional): output filename. Defaults to "structure.cif".
        """
        if True:
            _kkr2cif(self.structure, lattice=latticetype, fmt="cif",
                     outputfilename=outputfilename)
        else:
            _writer(self.structure, lattice=latticetype, fmt="cif",
                    outputfilename=outputfilename)

    def poscar(self, latticetype="direct", outputfilename="POSCAR"):
        """write as poscar

        Args:
            latticetype (str, optional): lattice type. Defaults to "primtitive".
            outputfilename (str, optional): output filename. Defaults to "POSCAR".
        """
        _writer(self.structure, lattice=latticetype, fmt="poscar",
                outputfilename=outputfilename)

    def xsf(self, latticetype="standard", outputfilename="structure.xsf"):
        """write as xsf

        Args:
            latticetype (str, optional): xsf filename. Defaults to "standard".
            outputfilename (str, optional): output filename. Defaults to "structure.xsf".
        """
        _writer_xsf(self.structure, lattice=latticetype,
                    outputfilename=outputfilename)
