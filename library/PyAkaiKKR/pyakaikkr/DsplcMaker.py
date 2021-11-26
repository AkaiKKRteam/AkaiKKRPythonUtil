# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from copy import deepcopy
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core import Structure, Lattice
import numpy as np
import math
from pymatgen.core.sites import PeriodicSite
from scipy import integrate
from typing import List, Tuple
from pymatgen.core.operations import SymmOp
from sympy import sympify, symbols, Symbol

from .AkaiKkr import AkaikkrJob
from .KkrPrimVec import PrmrConverter


def _make_wyckoff_position_as_sympy(ops: List[SymmOp], xyz: List[str],
                                    value: list = [
                                        ("x", 1), ("y", 1), ("z", 1)],
                                    frac: bool = True):
    """make wyckoff positions from ops

    xyz can be ["x","x",0] and so on.

    Args:
        ops ([pymatgen.core.operations.SymmOp]): point group operations
        xyz ([str]): ["x","y","z"]
        value ([str,float]): values of xyz
        frac (bool, optional): fractional coordinate or not. Defaults to True.
    """

    print("debug, _make_wyckoff_position_as_sympy starts")

    def _op2stroplist(ops: List[SymmOp], xyz: List[str]) -> List[str]:
        """convert ops to a list of string operations

        Args:
            ops (List[SymmOp]): [description]
            xyz (List[str]): [description]

        Returns:
            List[str]: [description]
        """
        _xyz = ["x", "y", "z"]
        _oplist = []
        for _iop, _op in enumerate(ops):
            _xyzop = _op.as_xyz_string()
            if _xyzop not in _oplist:
                for _a, _sub in zip(_xyz, xyz):
                    _xyzop = _xyzop.replace(_a, _sub)
                _oplist.append(_xyzop)
        return _oplist

    def _op2symoplist(oplist: List[str]) -> list:
        """convert a list of string operation to sympy instances

        Args:
            oplist (List[str]): [description]

        Returns:
            list: [description]
        """
        _x, _y, _z = symbols("x y z")
        _symoplist = []
        for _op in oplist:
            _op3 = _op.split(",")
            _oplist = []
            for _op1 in _op3:
                _expr = sympify(_op1)
                _oplist.append(_expr)
            _symoplist.append(_oplist)
        return _symoplist

    def _symop2strsymop(symoplist: list) -> List[str]:
        """convert a list of sympy instances to string operations

        Args:
            symoplist (list): [description]

        Returns:
            List[str]: [description]
        """
        _strsymoplist = []
        for _symop3 in symoplist:
            _symoplist = []
            for _symop1 in _symop3:
                _symoplist.append(str(_symop1))
            _strsymop = ",".join(_symoplist)
            if _strsymop not in _strsymoplist:
                _strsymoplist.append(_strsymop)
        return _strsymoplist

    def _strsymop2symop(strsymoplist: List[str]) -> list:
        """convert a list of string operations to a list of sympy instances

        Args:
            strsymoplist (List[str]): [description]

        Returns:
            list: [description]
        """
        symoplist = []
        for symop3 in strsymoplist:
            symop3 = symop3.split(",")
            _symoplist = []
            for _symop1 in symop3:
                expr = sympify(_symop1)
                _symoplist.append(expr)
            symoplist.append(_symoplist)
        return symoplist

    def _symop2dsplc(symoplist: list, value: List[Tuple[str, float]]) -> list:
        """make values from symoplist

        Args:
            symoplist (list): [description]
            value (List[str,float]): [description]

        Returns:
            List[float,float,float]: [description]
        """
        _valuelist = []
        for _symop3 in symoplist:
            _v3 = []
            for _symop1 in _symop3:
                _v = _symop1.subs(value)
                _v3.append(_v)
            _valuelist.append(_v3)
        return _valuelist

    if not frac:
        raise ValueError("frac must be True.")
    oplist = _op2stroplist(ops, xyz)
    symoplist = _op2symoplist(oplist)
    strsymoplist = _symop2strsymop(symoplist)
    symoplist = _strsymop2symop(strsymoplist)
    print("debug symop", symoplist)

    dsplc = _symop2dsplc(symoplist, value)
    print("debug dsplc", dsplc)
    return dsplc


def _frac2cartcoord(lattice: Lattice, fraclist: list, norm: float) -> list:
    """convert fractional coordinate to cartesian coordnate

    Args:
        lattice (np.ndarray): lattice vectors
        fraclist (list): a list of fractional coordinate
        norm (float): magnitude of output vector

    Returns:
        list: a list of cartesian coordinate
    """
    # lattice = np.array(lattice)

    cartlist = []
    for frac in fraclist:
        frac = np.array(frac)
        # v = lattice @ frac
        # v = np.dot(frac, lattice)
        v = lattice.get_cartesian_coords(frac)
        #vnorm = np.linalg.norm(v)
        vnorm = math.sqrt(np.sum(v**2))
        v = v * (norm/vnorm)
        cartlist.append(v.tolist())
    return cartlist


def _frac2PeriodicSite(lattice: Lattice, fraclist: list, norm: float,
                       dummy_specie="H") -> List[PeriodicSite]:
    """convert fractional coordinate to cartesian coordnate

    Args:
        lattice (np.ndarray): lattice vectors
        fraclist (list): a list of fractional coordinate
        norm (float): magnitude of output vector

    Returns:
        [PeriodicSite]: a list of PeriodicSite
    """

    cartlist = []
    for frac in fraclist:
        frac = np.array(frac)

        # v = lattice @ frac
        # v = np.dot(frac, lattice)
        v = lattice.get_cartesian_coords(frac)
        #vnorm = np.linalg.norm(v)
        vnorm = math.sqrt(np.sum(v**2))
        v = v * (norm/vnorm)
        psite = PeriodicSite(species=dummy_specie, coords=v, lattice=lattice,
                             coords_are_cartesian=True)
        print("debug, psite", psite)
        cartlist.append(psite)
    return cartlist


class DsplcMaker:
    def __init__(self, structure: Structure, kkrstruc: dict):
        """initialization

        Args:
            structure (Structure): pymatgen Sructure
            kkrstruc (dict): kkr structure parameter
        """
        self.structure = structure
        self.kkrstruc = kkrstruc
        print("debug kkrstruc", kkrstruc)

    def get_displc_as_periodicsite(self,
                                   xyz: list,
                                   cart_len=0.01,
                                   frac=True,
                                   show_detail=False):
        """get a list of displc as [PeriodicSite]

        For example, xyz and xyzvalue is given by, 
            xyz = ["x","x","0"]

        Now frac must be True.

        exception of SpacegroupAnalyzer(struc) should be catched... 

        Args:
            xyz ([[str]]) : xyz
            cart_len (float, optional): lenght of dx in cartesian. Defaults to 0.01.
            frac (bool, optional): direction is given by fractional or not. Defaults to True.
            show_detail (bool, optional): print detail or not. Defaults to True.

        Returns:
            [PeriodicSite]: rotated/operated PeriodicSites
        """

        if not frac:
            raise ValueError("frac must be True.")

        struc = self.structure
        spg = SpacegroupAnalyzer(struc)
        if show_detail:
            print("spginfo", spg.get_space_group_symbol(),
                  spg.get_space_group_number(),)

        if True:
            ops = spg.get_point_group_operations()
        else:
            ops = spg.get_space_group_operations()

        operatedsites = []
        for _xyz in xyz:
            fraclist = _make_wyckoff_position_as_sympy(
                ops, _xyz, frac=frac)

            prmrconverter = PrmrConverter()
            prmvec = prmrconverter.apply(self.kkrstruc["brvtyp"],
                                         self.kkrstruc["c/a"], self.kkrstruc["b/a"],
                                         self.kkrstruc["alpha"], self.kkrstruc["beta"],
                                         self.kkrstruc["gamma"])

            lattice = Lattice(prmvec)
            _operatedsites = _frac2PeriodicSite(lattice, fraclist, cart_len)

            if show_detail:
                for newosite in _operatedsites:
                    if show_detail:
                        np.set_printoptions(suppress=True)
                        print("debug operated site", "frac", newosite.frac_coords,
                              "cart", newosite.coords,)

            operatedsites.extend(_operatedsites)

        return operatedsites

    def get_displc_as_cartesian(self, cart_len: float,
                                xyz: list,
                                alen: float,
                                frac=True, show_detail=False):
        """get a list of displc as a list of [float,float,float]

           displc is scaled by alen

        Args:
            cart_len (float): lenght of dx in cartesian. 
            xyz ([float,float,float]): ["x", "y", "z"].
            xyzvalue ([(str,float)]) : values of xyz
            alen (float): lattice constant .
            frac (bool, optional): direction is given by fractional or not. Defaults to True.
            show_detail (bool, optional): print detail or not. Defaults to True.

        Returns:
            [[float,float,float]]: rotated/operated cartesian coordinates
        """
        operatedsites = self.get_displc_as_periodicsite(
            cart_len=cart_len,
            xyz=xyz,
            frac=frac,
            show_detail=show_detail)
        coords = []
        for sites in operatedsites:
            v = sites.coords/alen
            coords.append(v.tolist())

        return coords

    def get_avr_debye(self, m, debye_temp, temp):
        """get average displacement correspondinig to Debye temperature

        Args:
            m (float): mass
            debye_temp (float): Debye temperature
            temp (float): temperature

        Returns:
            float: average displacement
        """
        # mean square dispalcement by debye approximation
        # define parameters
        kB = 6.3336e-06            # Boltzmann constant, Ry/K
        me = 9.10938356            # electron mass, kg (10**(-31))
        u_to_kg = 1.660540         # conversion u to kg (10**(-27))
        #
        # convert atomic mass unit to atomic unit (i.e., me=0.5)
        m = m*u_to_kg/(2.0*me)*1.0e4
        #
        # define parameters for numerical integration
        mesh = 10000               # mesh for integration
        xinit = 1e-6               # initial point for integration
        #
        # debey function
        tau = debye_temp/temp
        t = np.linspace(xinit, tau, mesh)
        f = t/(np.exp(t)-1)
        D1 = integrate.trapz(f, t)/tau
        # D1=integrate.simps(f,t)/tau
        #
        # mean square displacement
        fac = 9.0/(m*kB*debye_temp)
        u2 = fac*(D1/tau)  # zero point energy is ignored
        # u2=fac*(D1/tau+0.25)
        u = np.sqrt(u2)
        return u


if __name__ == "__main__":
    import pprint
    import os
    parent_dir = "/home/kino/kino/kit"
    # path_dir = os.path.join(
    #    parent_dir, "AkaiKKRPythonUtil/tests/akaikkr_cnd/Co2MnSi")
    path_dir = os.path.join(
        parent_dir, "AkaiKKRPythonUtil/tests/akaikkr_cnd/SmCo5_oc")
    # path_dir = os.path.join(parent_dir,"AkaiKKRPythonUtil/tests/akaikkr_cnd/Co")
    # path_dir = os.path.join(
    #    parent_dir, "AkaiKKRPythonUtil/tests/akaikkr_cnd/Fe")
    # path_dir = os.path.join(
    #     parent_dir, "AkaiKKRPythonUtil/tests/akaikkr_cnd/Cu")

    specx = "/home/kino/kit/Akaikkrprogram.current.gfortran/akaikkr_cnd/specx"
    outgo = "out_go.log"
    cart_len = 0.01

    job = AkaikkrJob(path_dir)
    alen = job.get_lattice_constant(outgo)
    struc = job.make_pymatgenstructure(outgo)
    struc_param = job.get_struc_param(outgo)

    dsplcmaker = DsplcMaker(struc, struc_param)

    show_detail = False
    frac = True

    if True:
        # "ab" direction
        displc = dsplcmaker.get_displc_as_cartesian(
            cart_len=0.01, alen=alen,
            xyz=[["x", "0", "0"]],
            frac=frac,
            show_detail=show_detail)
        print("displc list")
        pp = pprint.PrettyPrinter(indent=1)
        pp.pprint(displc)
    if True:
        # "c" direction
        displc = dsplcmaker.get_displc_as_cartesian(
            cart_len=0.01, alen=alen,
            xyz=[["0", "0", "x"]],
            frac=frac,
            show_detail=show_detail)
        print("displc list")
        pp = pprint.PrettyPrinter(indent=1)
        pp.pprint(displc)
    if True:
        # bcc direction
        displc = dsplcmaker.get_displc_as_cartesian(
            cart_len=0.01, alen=alen,
            xyz=[["x", "x", "x"]],
            frac=frac,
            show_detail=show_detail)
        print("displc list")
        pp = pprint.PrettyPrinter(indent=1)
        pp.pprint(displc)
