# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import pandas as pd
from copy import deepcopy
import datetime
import numpy as np
import sys
import argparse
import os
from collections import OrderedDict
import json
import subprocess

from pyakaikkr import AkaikkrJob, JijPlotter
from pyakaikkr.CompareCifKkr import CompareCifKkr
from pyakaikkr.GoGo import *

from .exeutil import ExeUtil
from .resultutil import ResultUtil
from .OutputAnalyzer import *
from .loadmeta import make_meta

_REFERENCE_PATH_ = "reference"
_CHANGE_TYPE_IN_CIF_ = True


def get_kkr_struc_from_cif(ciffilepath: str, akaikkr_exe: str, displc: bool,
                           use_bravais=True, remove_temperaryfiles=True,
                           Vc: str = "Og",
                           directory: str = "temporary") -> dict:
    """get kkr structure parameter from cif file.

    If specx is akaikkr, then displc is set to False.
    If specx is akaikkr_cnd, then displc is set to True.

    If use_bravais is True, bcc,bcc,.., and a,c,b,alpha,beta,gaam are used.
    If use_bravias is False, aux and lattice vectors iare used.

    rmt, lmax, ..., displc are set to be default values.

    Args:
        ciffilepath (str): cif file path
        akaikkr_exe (str, optional): specx path. Defaults to str.
        displc (bool, optional): displc parameter is added.
        use_bravais (bool, optional): use bravias lattice. Defaults to True.
        remove_temperaryfiles (bool, optional): delete temporary files on exit. Defaults to True.
        Vc (str, optional): dummy vacancy element. Defaults to "Og"
        parent_directory (str, optional): working directory. Defaults to "temporary".

    Returns:
        dict: kkr structure parameters on success.
    """
    comp = CompareCifKkr(ciffilepath, akaikkr_exe, displc=displc,
                         Vc=Vc, parent_directory=directory)
    result = comp.convert_and_compare(use_bravais=use_bravais)
    struc_param = None
    if result == comp.SUCCESS:
        struc_param = comp.get_structure_param(
            remove_temperaryfiles=remove_temperaryfiles)
    else:
        print("failed to convert the cif file")
        print("msg=", comp.msg)
        print("result=", result)
    try:
        os.rmdir(comp.parent_directory)
    except OSError:
        # ignore errors on removing output directory
        pass
    return struc_param


def _atmicx_float2frac(atmicx: list, eps=1e-7) -> list:
    """convert 0.333333 to 1/3, also for 2/3

    Args:
        atmicx (list): atmicx
        eps (float, optional): eps to compare numerical values. Defaults to 1e-7.

    Returns:
        list: atmicx with fractional coordiates when they are 1/2 or 2/3
    """
    newatomicx = []
    for a_atmicx in atmicx:
        new_a_atmicx = []
        for p in a_atmicx[:3]:
            if p[-1] == "a" or p[-1] == "b" or p[-1] == "c":
                val = float(p[:-1])
                axis = p[-1]
            else:
                val = float(p)
                axis = ""
            newvalue = None
            for value, strvalue in zip([0.33333333333, 0.6666666666666], ["1/3", "2/3"]):
                if np.abs(val-value) < eps:
                    newvalue = "{}{}".format(strvalue, axis)
            if newvalue is None:
                newvalue = "{}{}".format(val, axis)
            new_a_atmicx.append(newvalue)
        new_a_atmicx.append(a_atmicx[-1])
        newatomicx.append(new_a_atmicx)
    return newatomicx


def _sep_atmicx(atmicx1: list) -> list:
    """convert atmicx format

    0.1a 0.2b 0.3c Mn -> 0.1 a 0.2 b 0.3 c Mn
    or
    0.1 0.2 0.3 Mn -> 0.1 0.2 0.3 Mn

    Args:
        atmicx1 (list): list of atmicx

    Returns:
        list: list of atmicx
    """
    newatmicx1 = []
    for a_atmicx in atmicx1:
        a_vector = []
        for p in a_atmicx[:3]:
            if p[-1] == "a" or p[-1] == "b" or p[-1] == "c":
                try:
                    v = float(p[:-1])
                except ValueError:
                    v = p[:-1]
                a_vector.extend([v, p[-1]])
            else:
                a_vector.append(p)
        a_vector.append(a_atmicx[-1])
        newatmicx1.append(np.array(a_vector))
    return newatmicx1


def _reorder_atmicx(atmicx1: np.ndarray, atmicx2: np.ndarray) -> np.ndarray:
    """reorder atmicx2 as the order of atmicx1

    Args:
        atmicx1 (list): atmicx1
        atmicx2 (list): atmicx2

    Raises:
        ValueError: len(newatmicx) != len(atmicx1)

    Returns:
        list: reordered atmicx2 as the order of atmicx1
    """
    newatmicx1 = _sep_atmicx(atmicx1)
    newatmicx2 = _sep_atmicx(atmicx2)
    newatmicx = []
    for atm1 in newatmicx1:
        for atm2 in newatmicx1:
            if np.all(atm1 == atm2):
                newatmicx.append(atm2.tolist())
                break
    if len(newatmicx) != len(atmicx1):
        print(newatmicx1)
        print(newatmicx2)
        print(newatmicx)
        raise ValueError("len(newatmicx) != len(atmicx1)")

    newatmicx_fmt = []
    for atm in newatmicx:
        newatmicx_fmt.append(["{}{}".format(atm[0], atm[1]),
                              "{}{}".format(atm[2], atm[3]),
                              "{}{}".format(atm[4], atm[5]),
                              atm[6]])
    return newatmicx_fmt


def _reorder_type(param: dict, typeorder: list) -> dict:
    """reorder param by typeorder

    Args:
        param (dict): kkr pamarameter
        typeorder ([str]): a list of types

    Returns:
        dict: reordered kkr parameter
    """
    rmt = []
    field = []
    conc = []
    anclr = []
    mxl = []
    ncmp = []
    for newtype in typeorder:
        i = param["type"].index(newtype)
        rmt.append(param["rmt"][i])
        field.append(param["field"][i])
        conc.append(param["conc"][i])
        anclr.append(param["anclr"][i])
        mxl.append(param["mxl"][i])
        ncmp.append(param["ncmp"][i])

    param["type"] = typeorder
    param["rmt"] = rmt
    param["field"] = field
    param["conc"] = conc
    param["anclr"] = anclr
    param["mxl"] = mxl
    param["ncmp"] = ncmp
    return param


def change_atomic_type(param: dict, type_rep: dict) -> dict:
    """change type by type_rep.

    change type param["type"] and param["atmicx"] according to type_rep = {"Cu_4a_0": "Cu"}.

    Args:
        param (dict): kkr param
        type_rep (dict): replace dict
    """
    newtype = []
    for type_ in param["type"]:
        newtype.append(type_rep[type_])
    param["type"] = newtype
    newatmicx = []
    for atmicx in param["atmicx"]:
        newatmicx.append(
            [atmicx[0], atmicx[1], atmicx[2], type_rep[atmicx[3]]])
    param["atmicx"] = newatmicx
    return param


def _Cu_common_param(
        akaikkr_exe: dict, displc: bool, ciffilepath="../structure/Cu-Fm3m.cif",
        use_bravais=True, remove_temperaryfiles=False,
        directory="temporary") -> dict:
    """define Cu common kkr parameters

    Args:
        akaikkr_exe (dict): specx and fmg
        displc (bool): include displc or not. True if akaikkr_cnd.
        ciffilepath (str, optional): cif file name. Defaults to "../structure/Cu-Fm3m.cif".
        use_bravais (bool, optional): use bravias in the akaikkr format. Defaults to True.
        remove_temperaryfiles (bool, optional): remove temporary files after executing geometry by specx. Defaults to False.
        directory (str, optional): directory to place input and output files in getting geometry. Defaults to "temporary".

    Returns:
        dict: kkr input parameters
    """

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 6.82
    param["magtyp"] = "nmag"

    # to compare the results with reference
    param = change_atomic_type(param, {"Cu_4a_0": "Cu"})

    return param


_Cu_COMMENT_ = "fcc Cu (non-magnetic, mjw, {})"


def Cu_go(akaikkr_exe: dict, directory="Cu", comment=_Cu_COMMENT_,
          displc=False, execute_postscript=True):
    """define Cu go execution

    Args:
        akaikkr_exe (dict): define specx and fmg
        directory (str, optional): directory to place input and output of akaikkr. Defaults to "Cu".
        comment (str, optional): short description of this test. Defaults to _Cu_COMMENT_.
        displc (bool, optional): include displc or not. Defaults to False.
        execute_postscript (bool, optional): also execute visualization. Defaults to True.

    Returns:
        str, dict: label, result in the dict format
    """
    # option = {"ie":1,"mse":3}
    # args = {"option": option}
    gogo = GoGo(akaikkr_exe, directory,
                _Cu_common_param(akaikkr_exe=akaikkr_exe,
                                 displc=displc, directory=directory),
                comment)  # , **args)

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Cu_dos(akaikkr_exe, directory="Cu", comment=_Cu_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _Cu_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Cu_spc(akaikkr_exe,  directory="Cu", comment=_Cu_COMMENT_,
           displc=False, execute_postscript=True):
    # options = {"ie":1,"mse":3, "klabel": ["A","B","C"]}
    # args = {"options": options}
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _Cu_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, go=go)  # , **args)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _Fe_common_param(
        akaikkr_exe: dict, displc: bool, ciffilepath="../structure/Fe-Im3m.cif",
        use_bravais=True, remove_temperaryfiles=False,
        directory="temporary") -> dict:

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 5.27
    param["magtyp"] = "mag"
    param["rmt"] = [1.0 for i in param["ncmp"]]

    return param


_Fe_COMMENT_ = "bcc Fe (magnetic, mjw, {})"


def Fe_go(akaikkr_exe,  directory="Fe", comment=_Fe_COMMENT_,
          displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory, _Fe_common_param(akaikkr_exe=akaikkr_exe,
                                                         displc=displc, directory=directory),  comment, )

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_fsm(akaikkr_exe, directory="Fe", comment=_Fe_COMMENT_,
           displc=False, execute_postscript=True):
    _param_add = {"potentialfile": "pot_fsm.dat", "record": "init"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _Fe_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_tc(akaikkr_exe,  directory="Fe", comment=_Fe_COMMENT_,
          displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,
                _Fe_common_param(akaikkr_exe=akaikkr_exe,
                                 displc=displc, directory=directory), comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_j30(akaikkr_exe, directory="Fe", comment=_Fe_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,
                 _Fe_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory), comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    if execute_postscript:
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_typepair()
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_dos(akaikkr_exe, directory="Fe", comment=_Fe_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _Fe_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory), comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_spc(akaikkr_exe, directory="Fe",  comment=_Fe_COMMENT_,
           displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _Fe_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory), comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _Co_common_param(
    akaikkr_exe: dict, displc: bool, ciffilepath="../structure/Co_P63mmc.cif",
        use_bravais=True, remove_temperaryfiles=False,
        directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 4.74
    # param["c/a"] = 1.6215
    param["magtyp"] = "mag"
    param["rmt"] = [1.0 for i in param["ncmp"]]

    return param


_Co_COMMENT_ = "hcp Co (magnetic, mjw, {})"


def Co_go(akaikkr_exe, directory="Co", comment=_Co_COMMENT_,
          displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _Co_common_param(akaikkr_exe=akaikkr_exe,
                                 displc=displc, directory=directory),
                comment, )
    gogo.execute(displc=displc, execute_postscript=True)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co_fsm(akaikkr_exe,  directory="Co", comment=_Co_COMMENT_,
           displc=False, execute_postscript=True):
    _param_add = {"potentialfile": "pot_fsm.dat", "record": "init"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _Co_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co_tc(akaikkr_exe, directory="Co", comment=_Co_COMMENT_,
          displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,
                _Co_common_param(akaikkr_exe=akaikkr_exe,
                                 displc=displc, directory=directory),
                comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co_j30(akaikkr_exe,  directory="Co", comment=_Co_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,
                 _Co_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    if execute_postscript:
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_typepair()
    return label, gogo.result


def Co_dos(akaikkr_exe, directory="Co", comment=_Co_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _Co_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co_spc(akaikkr_exe, directory="Co", comment=_Co_COMMENT_,
           displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _Co_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def Co_spc21(akaikkr_exe, directory="Co", comment=_Co_COMMENT_,
             displc=False, execute_postscript=True):
    go = "spc21"
    args = {"fmt": 2, "first_connected_kpath": False}
    gogo = GoSpc(akaikkr_exe, directory,
                 _Co_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, go=go, **args)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _Ni_common_param(akaikkr_exe: dict, displc: bool, ciffilepath="../structure/Ni-Fm3m.cif",
                     use_bravais=True, remove_temperaryfiles=False,
                     directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    param["magtyp"] = "mag"
    param["rmt"] = [1.0 for i in param["ncmp"]]

    return param


_Ni_COMMENT_ = "fcc Ni (magnetic, mjw, {})"


def Ni_go(akaikkr_exe,  directory="Ni", comment=_Ni_COMMENT_,
          displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _Ni_common_param(akaikkr_exe=akaikkr_exe,
                                 displc=displc, directory=directory),
                comment, )

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Ni_fsm(akaikkr_exe,  directory="Ni", comment=_Ni_COMMENT_,
           displc=False, execute_postscript=True):
    _param_add = {"potentialfile": "pot_fsm.dat", "record": "init"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _Ni_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Ni_tc(akaikkr_exe,  directory="Ni", comment=_Ni_COMMENT_,
          displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,
                _Ni_common_param(akaikkr_exe=akaikkr_exe,
                                 displc=displc, directory=directory),
                comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Ni_j30(akaikkr_exe,  directory="Ni", comment=_Ni_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,
                 _Ni_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    if execute_postscript:
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_typepair()
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Ni_dos(akaikkr_exe,  directory="Ni", comment=_Ni_COMMENT_,
           displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _Ni_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Ni_spc(akaikkr_exe, directory="Ni", comment=_Ni_COMMENT_,
           displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _Ni_common_param(akaikkr_exe=akaikkr_exe,
                                  displc=displc, directory=directory),
                 comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _AlMnFeCo_bcc_common_param(akaikkr_exe: dict, displc: bool, ciffilepath="../structure/AlMnFeCo-Im3m.cif",
                               use_bravais=True, remove_temperaryfiles=False,
                               directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    param["a"] = 1000000
    param["magtyp"] = "mag"
    param["reltyp"] = "sra"
    param["magtyp"] = "mag"
    param["sdftyp"] = "pbeasa"
    param["ewidth"] = 1.3
    param["pmix"] = 0.01
    param["rmt"] = [1.0]
    param["mxl"] = [3]

    return param


_AlMnFeCo_bcc_COMMENT_ = "bcc AlMnFeCo (magnetic, pbeasa, {})"


def AlMnFeCo_bcc_go(akaikkr_exe,  directory="AlMnFeCo_bcc", comment=_AlMnFeCo_bcc_COMMENT_,
                    displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                           displc=displc, directory=directory),
                comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def AlMnFeCo_bcc_gofmg(akaikkr_exe,  directory="AlMnFeCo_bcc",
                       comment="bcc AlMnFeCo (magnetic, pbeasa, spin flipped Mn, {})",
                       displc=False, execute_postscript=True):

    flip_list = ["HEA_Mn_25.0%", ]
    gogo = GoFmg(akaikkr_exe, directory,
                 _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                            displc=displc, directory=directory),
                 comment,
                 flip_list,)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, "gofmg")
    return label, gogo.result


def AlMnFeCo_bcc_fsm(akaikkr_exe, directory="AlMnFeCo_bcc", comment=_AlMnFeCo_bcc_COMMENT_,
                     displc=False, execute_postscript=True):
    _param_add = {"potentialfile": "pot_fsm.dat", "record": "init"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                            displc=displc, directory=directory),
                 comment, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def AlMnFeCo_bcc_tc(akaikkr_exe, directory="AlMnFeCo_bcc", comment=_Ni_COMMENT_,
                    displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,
                _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                           displc=displc, directory=directory),
                comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def AlMnFeCo_bcc_j30(akaikkr_exe,  directory="AlMnFeCo_bcc", comment=_AlMnFeCo_bcc_COMMENT_,
                     displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,
                 _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                            displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    if execute_postscript:
        job = AkaikkrJob(directory)
        outfile = "out_go.log"
        typeofsite = job.get_type_of_site(outfile)
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_comppair(
            "Mn0.25Al0.25Fe0.25Co0.25_2a_0", "Mn0.25Al0.25Fe0.25Co0.25_2a_0", typeofsite,
        )
    return label, gogo.result


def AlMnFeCo_bcc_dos(akaikkr_exe,  directory="AlMnFeCo_bcc", comment=_AlMnFeCo_bcc_COMMENT_,
                     displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                            displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def AlMnFeCo_bcc_cnd(akaikkr_exe, directory="AlMnFeCo_bcc", comment=_AlMnFeCo_bcc_COMMENT_,
                     displc=False, execute_postscript=True):
    gogo = GoCnd(akaikkr_exe, directory,
                 _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                            displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def AlMnFeCo_bcc_spc(akaikkr_exe, directory="AlMnFeCo_bcc", comment=_AlMnFeCo_bcc_COMMENT_,
                     displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _AlMnFeCo_bcc_common_param(akaikkr_exe=akaikkr_exe,
                                            displc=displc, directory=directory),
                 comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _FeRh05Pt05_common_param(akaikkr_exe: dict, displc: bool, ciffilepath="../structure/FeRh0.5Pt0.5.cif",
                             use_bravais=True, remove_temperaryfiles=False,
                             directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 5.2043
    # param["c/a"] = 1.2817
    param["magtyp"] = "mag"
    param["sdftyp"] = "pbeasa"
    param["ewidth"] = 1.0
    param["pmix"] = "0.02ch"
    param["bzqlty"] = 8
    param["rmt"] = [1.0 for i in range(param["ntyp"])]
    param["mxl"] = [3 for i in range(param["ntyp"])]

    return param


_FeRh05Pt05_COMMENT_ = "so FeRh0.5Pt0.5. (magnetic, pbeasa, {})"


def FeRh05Pt05_go(akaikkr_exe, directory="FeRh05Pt05",
                  comment=_FeRh05Pt05_COMMENT_,
                  displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                         displc=displc, directory=directory),
                comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeRh05Pt05_fsm(akaikkr_exe,  directory="FeRh05Pt05",
                   comment=_FeRh05Pt05_COMMENT_,
                   displc=False, execute_postscript=True):
    # copy potential and change fspin a little
    job = AkaikkrJob(directory)
    file1 = job.default["potentialfile"]
    file2 = "pot_fsm.dat"
    job.copy_potential(file1, file2)
    fspin = 3.0
    _param_add = {"potentialfile": file2, "record": "2nd"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                          displc=displc, directory=directory),
                 comment, fspin=fspin, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeRh05Pt05_dos(akaikkr_exe, directory="FeRh05Pt05",
                   comment=_FeRh05Pt05_COMMENT_,
                   displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory, _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeRh05Pt05_tc(akaikkr_exe,  directory="FeRh05Pt05",
                  comment=_FeRh05Pt05_COMMENT_,
                  displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory, _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                                                 displc=displc, directory=directory),
                comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeRh05Pt05_j30(akaikkr_exe,  directory="FeRh05Pt05",
                   comment=_FeRh05Pt05_COMMENT_,
                   displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory, _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    if execute_postscript:
        outfile = "out_go.log"
        job = AkaikkrJob(directory)
        typeofsite = job.get_type_of_site(outfile)
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_comppair(
            "Fe_1a_0", "Rh0.5Pt0.5_1d_1", typeofsite, )
        jijplotter.make_comppair("Fe_1a_0", "Fe_1a_0",
                                 typeofsite,)

    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeRh05Pt05_cnd(akaikkr_exe, directory="FeRh05Pt05",
                   comment=_FeRh05Pt05_COMMENT_,
                   displc=False, execute_postscript=True):
    gogo = GoCnd(akaikkr_exe, directory, _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeRh05Pt05_spc(akaikkr_exe,  directory="FeRh05Pt05",
                   comment=_FeRh05Pt05_COMMENT_,
                   displc=False, execute_postscript=True):
    gogo = GoSpc(akaikkr_exe, directory, _FeRh05Pt05_common_param(akaikkr_exe=akaikkr_exe,
                                                                  displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def _NiFe_common_param(akaikkr_exe: dict, displc: bool, ciffilepath="../structure/NiFe-Fm3m.cif",
                       use_bravais=True, remove_temperaryfiles=False,
                       directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 6.65
    param["magtyp"] = "mag"
    param["sdftyp"] = "pbeasa"
    param["bzqlty"] = 8  # =6 is NG
    if param["ntyp"] != 1:
        raise ValueError("ntyp must be 1 for NiFe.")
    param["rmt"] = [1.0 for i in range(param["ntyp"])]
    param["mxl"] = [3 for i in range(param["ntyp"])]

    return param


_NiFe_COMMENT_ = "fcc NiFe (magnetic, pbeasa, {})"


def NiFe_go(akaikkr_exe, directory="NiFe", comment=_NiFe_COMMENT_,
            displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                   displc=displc, directory=directory),
                comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def NiFe_fsm(akaikkr_exe,  directory="NiFe", comment=_NiFe_COMMENT_,
             displc=False, execute_postscript=True):
    _param_add = {"potentialfile": "pot_fsm.dat", "record": "init"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def NiFe_tc(akaikkr_exe,  directory="NiFe", comment=_NiFe_COMMENT_,
            displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,
                _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                   displc=displc, directory=directory),
                comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def NiFe_j30(akaikkr_exe,  directory="NiFe", comment=_NiFe_COMMENT_,
             displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,
                 _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    if execute_postscript:
        outfile = "out_go.log"
        job = AkaikkrJob(directory)
        typeofsite = job.get_type_of_site(outfile)
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_comppair(
            "Fe0.1Ni0.9_4a_0", "Fe0.1Ni0.9_4a_0", typeofsite, )
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def NiFe_dos(akaikkr_exe, directory="NiFe", comment=_NiFe_COMMENT_,
             displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def NiFe_cnd(akaikkr_exe,  directory="NiFe", comment=_NiFe_COMMENT_,
             displc=False, execute_postscript=True):
    gogo = GoCnd(akaikkr_exe, directory,
                 _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def NiFe_spc(akaikkr_exe, directory="NiFe", comment=_NiFe_COMMENT_,
             displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _NiFe_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _Fe_lmd_common_param(akaikkr_exe: dict, displc: bool, ciffilepath="../structure/Fe-Im3m.cif",
                         use_bravais=True, remove_temperaryfiles=False,
                         directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 5.27
    param["magtyp"] = "lmd"
    param["rmt"] = [1.0 for i in range(param["ntyp"])]

    # lmd must increase ncomp
    ncmp_list = []
    for ncmp in param["ncmp"]:
        if not isinstance(ncmp, int):
            raise ValueError("Each ncmp must be type int.")
        ncmp_list.append([ncmp, ncmp])
    param["ncmp"] = [ncmp*2]

    # alloy isn't allowed now.
    # extend sizes twice
    anclr_list = []
    for anclr in param["anclr"]:
        _x = anclr
        if len(_x) != 1:
            raise ValueError("Each anclr must be 1.")
        _x.extend(_x)
        anclr_list.append(_x)
    param["anclr"] = anclr_list

    conc_list = []
    for conc in param["conc"]:
        conc = np.array(conc)
        _x = conc*0.5
        _x = _x.tolist()
        if len(_x) != 1:
            raise ValueError("Each conc must be 1.")
        _x.extend(_x)
        conc_list.append(_x)
    param["conc"] = conc_list

    return param


_Fe_lmd_COMMENT_ = "bcc Fe(lmd, mjw, {})"


def Fe_lmd_go(akaikkr_exe, directory="Fe_lmd", comment=_Fe_lmd_COMMENT_,
              displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory, _Fe_lmd_common_param(akaikkr_exe=akaikkr_exe,
                                                             displc=displc, directory=directory),  comment, )

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_lmd_tc(akaikkr_exe, directory="Fe_lmd", comment=_Fe_lmd_COMMENT_,
              displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,  _Fe_lmd_common_param(akaikkr_exe=akaikkr_exe,
                                                              displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_lmd_j30(akaikkr_exe,  directory="Fe_lmd", comment=_Fe_lmd_COMMENT_,
               displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,  _Fe_lmd_common_param(akaikkr_exe=akaikkr_exe,
                                                               displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_lmd_dos(akaikkr_exe,  directory="Fe_lmd", comment=_Fe_lmd_COMMENT_,
               displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,  _Fe_lmd_common_param(akaikkr_exe=akaikkr_exe,
                                                               displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Fe_lmd_spc(akaikkr_exe,  directory="Fe_lmd", comment=_Fe_lmd_COMMENT_,
               displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _Fe_lmd_common_param(akaikkr_exe=akaikkr_exe,
                                      displc=displc, directory=directory),  comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _FeB195_common_param(akaikkr_exe: dict, displc: bool,
                         ciffilepath="../structure/FeB1.95-P6mmm.cif",
                         use_bravais=True, remove_temperaryfiles=False,
                         directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        Vc="Og", directory=directory)

    param["edelt"] = 1e-3
    param["ewidth"] = 1.5
    param["magtyp"] = "mag"

    return param


_FeB195_COMMENT_ = "FeB1.95 (magnetic, mjw, {})"


def FeB195_go(akaikkr_exe,  directory="FeB195", comment=_FeB195_COMMENT_,
              displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _FeB195_common_param(akaikkr_exe=akaikkr_exe,
                                     displc=displc, directory=directory),
                comment,)

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeB195_dos(akaikkr_exe,  directory="FeB195", comment=_FeB195_COMMENT_,
               displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _FeB195_common_param(akaikkr_exe=akaikkr_exe,
                                      displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def FeB195_spc(akaikkr_exe,  directory="FeB195", comment=_FeB195_COMMENT_,
               displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _FeB195_common_param(akaikkr_exe=akaikkr_exe,
                                      displc=displc, directory=directory),  comment, go=go)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _GaAs_common_param(akaikkr_exe: dict, displc: bool,
                       ciffilepath="../structure/GaAsVc-F43m.cif",
                       use_bravais=True, remove_temperaryfiles=False,
                       directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        Vc="Og", directory=directory)

    # param["a"] = 10.684
    param["edelt"] = 1e-3
    param["ewidth"] = 1.5
    param["magtyp"] = "nmag"

    param["ncmp"] = []
    param["rmt"] = []
    param["mxl"] = []
    for ityp in range(param["ntyp"]):
        param["ncmp"].append(1)
        param["rmt"].append(1.0)
        param["mxl"].append(2)

    return param


_GaAs_COMMENT_ = "GaAs (non-magnetic, mjw, {})"


def GaAs_go(akaikkr_exe,  directory="GaAs", comment=_GaAs_COMMENT_,
            displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _GaAs_common_param(akaikkr_exe=akaikkr_exe,
                                   displc=displc, directory=directory),
                comment,)

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def GaAs_dos(akaikkr_exe, directory="GaAs", comment=_GaAs_COMMENT_,
             displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,
                 _GaAs_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment,)

    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def GaAs_spc(akaikkr_exe,  directory="GaAs", comment=_GaAs_COMMENT_,
             displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory,
                 _GaAs_common_param(akaikkr_exe=akaikkr_exe,
                                    displc=displc, directory=directory),
                 comment,)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _Co2MnSi_common_param(akaikkr_exe: dict, displc: bool,
                          ciffilepath="../structure/Co2MnSi-Fm3m.cif",
                          use_bravais=True, remove_temperaryfiles=False,
                          directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 10.6675
    param["edelt"] = 1e-4
    param["ewidth"] = 1.2
    param["magtyp"] = "mag"
    param["sdftyp"] = "pbe"
    param["rmt"] = [1.0 for i in range(param["ntyp"])]

    return param


_Co2MnSi_COMMENT_ = "Co2MnSi (magnetic, pbe, {})"


def Co2MnSi_go(akaikkr_exe,  directory="Co2MnSi", comment=_Co2MnSi_COMMENT_,
               displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _Co2MnSi_common_param(akaikkr_exe=akaikkr_exe,
                                      displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co2MnSi_fsm(akaikkr_exe,  directory="Co2MnSi", comment=_Co2MnSi_COMMENT_,
                displc=False, execute_postscript=True):
    # copy potential and change fspin a little
    job = AkaikkrJob(directory)
    file1 = job.default["potentialfile"]
    file2 = "pot_fsm.dat"
    job.copy_potential(file1, file2)
    fspin = 4.5
    _param_add = {"potentialfile": file2, "record": "2nd"}
    gogo = GoFsm(akaikkr_exe, directory,
                 _Co2MnSi_common_param(akaikkr_exe=akaikkr_exe,
                                       displc=displc, directory=directory),  comment, fspin=fspin, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co2MnSi_tc(akaikkr_exe,  directory="Co2MnSi", comment=_Co2MnSi_COMMENT_,
               displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,  _Co2MnSi_common_param(akaikkr_exe=akaikkr_exe,
                                                               displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co2MnSi_j30(akaikkr_exe,  directory="Co2MnSi", comment=_Co2MnSi_COMMENT_,
                displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,  _Co2MnSi_common_param(akaikkr_exe=akaikkr_exe,
                                                                displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    if execute_postscript:
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_typepair()
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co2MnSi_dos(akaikkr_exe,  directory="Co2MnSi", comment=_Co2MnSi_COMMENT_,
                displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,  _Co2MnSi_common_param(akaikkr_exe=akaikkr_exe,
                                                                displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def Co2MnSi_spc(akaikkr_exe,  directory="Co2MnSi", comment=_Co2MnSi_COMMENT_,
                displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory, _Co2MnSi_common_param(akaikkr_exe=akaikkr_exe,
                                                               displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)

    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _SmCo5_oc_common_param(akaikkr_exe: dict, displc: bool, ciffilepath="../structure/SmCo5_P6mmm.cif",
                           use_bravais=True, remove_temperaryfiles=False,
                           directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    # param["a"] = 9.415
    # param["c/a"] = 0.798

    param["ewidth"] = 1.7
    param["reltyp"] = "srals"
    param["sdftyp"] = "mjwasa"
    param["magtyp"] = "mag"
    param["rmt"] = [0.0 for i in range(param["ntyp"])]

    return param


_SmCo5_oc_COMMENT_ = "SmCo5(magnetic, open core for Sm-4f, srals, mjw-asa, {})"


def SmCo5_oc_go(akaikkr_exe,  directory="SmCo5_oc", comment=_SmCo5_oc_COMMENT_,
                displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory, _SmCo5_oc_common_param(akaikkr_exe=akaikkr_exe,
                                                               displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def SmCo5_oc_fsm(akaikkr_exe, directory="SmCo5_oc", comment=_SmCo5_oc_COMMENT_,
                 displc=False, execute_postscript=True):
    # copy potential and change fspin a little
    job = AkaikkrJob(directory)
    file1 = job.default["potentialfile"]
    file2 = "pot_fsm.dat"
    job.copy_potential(file1, file2)
    _param_add = {"potentialfile": file2, "record": "2nd"}
    fspin = 6.5
    gogo = GoFsm(akaikkr_exe, directory,
                 _SmCo5_oc_common_param(akaikkr_exe=akaikkr_exe,
                                        displc=displc, directory=directory),  comment, fspin=fspin, **_param_add)
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def SmCo5_oc_tc(akaikkr_exe,  directory="SmCo5_oc", comment=_SmCo5_oc_COMMENT_,
                displc=False, execute_postscript=True):
    gogo = GoTc(akaikkr_exe, directory,  _SmCo5_oc_common_param(akaikkr_exe=akaikkr_exe,
                                                                displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def SmCo5_oc_j30(akaikkr_exe,  directory="SmCo5_oc", comment=_SmCo5_oc_COMMENT_,
                 displc=False, execute_postscript=True):
    gogo = Goj30(akaikkr_exe, directory,  _SmCo5_oc_common_param(akaikkr_exe=akaikkr_exe,
                                                                 displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    if execute_postscript:
        jijplotter = JijPlotter(gogo.df_jij, directory)
        jijplotter.make_typepair()
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def SmCo5_oc_dos(akaikkr_exe,   directory="SmCo5_oc", comment=_SmCo5_oc_COMMENT_,
                 displc=False, execute_postscript=True):
    gogo = GoDos(akaikkr_exe, directory,  _SmCo5_oc_common_param(akaikkr_exe=akaikkr_exe,
                                                                 displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def SmCo5_oc_spc(akaikkr_exe,  directory="SmCo5_oc", comment=_SmCo5_oc_COMMENT_,
                 displc=False, execute_postscript=True):
    go = "spc31"
    gogo = GoSpc(akaikkr_exe, directory, _SmCo5_oc_common_param(akaikkr_exe=akaikkr_exe,
                                                                displc=displc, directory=directory),  comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, go)
    return label, gogo.result


def _SmCo5_noc_common_param(akaikkr_exe: dict, displc: bool,
                            ciffilepath="../structure/SmCo5_P6mmm.cif",
                            use_bravais=True, remove_temperaryfiles=False,
                            directory="temporary"):

    param = get_kkr_struc_from_cif(
        ciffilepath=ciffilepath, akaikkr_exe=akaikkr_exe["specx"], displc=displc,
        use_bravais=use_bravais, remove_temperaryfiles=remove_temperaryfiles,
        directory=directory)

    #param["a"] = 9.415
    #param["c/a"] = 0.798

    param["ewidth"] = "1.7noc"
    param["reltyp"] = "srals"
    param["sdftyp"] = "mjwasa"
    param["magtyp"] = "mag"
    param["rmt"] = [0.0 for i in range(param["ntyp"])]

    return param


_SmCo5_noc_COMMENT_ = "SmCo5 (magnetic, srals, mjw-asa, {})"


def SmCo5_noc_go(akaikkr_exe,  directory="SmCo5_noc", comment=_SmCo5_noc_COMMENT_,
                 displc=False, execute_postscript=True):
    gogo = GoGo(akaikkr_exe, directory,
                _SmCo5_noc_common_param(akaikkr_exe=akaikkr_exe,
                                        displc=displc, directory=directory),
                comment, )
    gogo.execute(displc=displc, execute_postscript=execute_postscript)
    label = "{}_{}".format(directory, gogo.go)
    return label, gogo.result


def parse_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        "program_path", help="path to specx and fmg. the path is without the last akaikkr or akaikkr_cnd")
    parser.add_argument(
        "--create_ref", help="create reference data", action="store_true")
    parser.add_argument(
        "--compiler", help="compiler of the reference data", default="ifort")
    parser.add_argument(
        "--no_run", help="not run, only read output", action="store_true")
    parser.add_argument(
        "--no_postscript", help="don't make images", action="store_true")
    parser.add_argument(
        "--set", help="specify exe_list set", default="all")
    args = parser.parse_args()
    return args


def get_commit_id():
    s = subprocess.check_output(["git", "log"])
    lines = s.splitlines()
    for line in lines:
        lined = line.decode()
        s = lined.split()
        if s[0] == "commit":
            commit = s[1]
            break
    return {"commit": commit}


def check_go(key: str, result_dic: dict, ref_dic_result: dict):
    go = result_dic[key]["go"]
    if ref_dic_result is None:
        ref = None
    else:
        if key in ref_dic_result:
            ref = ref_dic_result[key]
        else:
            ref = None

    if go == "go" or go == "fsm":
        chk = go_diff_msg(
            key, result_dic[key], ref)
    elif go == "tc":
        chk = tc_diff_msg(
            key, result_dic[key], ref)
    elif go[0] == "j":
        chk = j_diff_msg(
            key, result_dic[key], ref)
    elif go == "dos":  # totalDOS
        chk = dos_diff_msg(
            key, result_dic[key], ref)
    elif go == "cnd":
        chk = cnd_diff_msg(
            key, result_dic[key], ref)
    elif go[:3] == "spc":
        chk = aw_diff_msg(
            key, result_dic[key], ref)

    else:
        print("unknown go", go)
        raise ValueError("unknown go")
    return chk


def all_go(akaikkr_exe, fmg_exe, exe_dic, displc=False):
    # stop using it because line breakig occurs...
    # pd.set_option("display.max_columns", 12)
    pd.set_option("display.precision", 8)
    pd.options.display.float_format = '{:.2e}'.format

    args = parse_args()

    ref_file = os.path.join(_REFERENCE_PATH_, args.compiler+".json")
    print(ref_file)
    meta_file = "meta.json"
    meta_str = """{
"comment": "Reference calculations are performed on iceberg in NIMS.",
"author": "Hiori Kino",
"CPU": "Intel(R) Xeon(R) CPU E5-2640 v4 @ 2.40GHz",
"compiler": "ifort v16.0.2",
"option": "-O2 -mcmodel=medium -qopenmp"
}
"commit" and "date" fields are added automatically.
"""
    if not os.path.isfile(ref_file) and not args.create_ref:
        print("No reference data found.")
        print("Add --create_ref and run")
        print(meta_file, "must be prepared if you create reference data.")
        print(meta_file, "can be made automatically.")
        print("You can prepare meta.json by hand. An example of the content is")
        print(meta_str)
        sys.exit(1)
    if os.path.isfile(ref_file) and args.create_ref:
        print(ref_file, "exists. delete it first for --create_ref")
        sys.exit(1)

    if os.path.isfile(meta_file):
        print("reading", meta_file)
        with open(meta_file) as f:
            meta = json.load(f)
            now = datetime.datetime.utcnow().ctime()
            print("datetime is updated.")
            meta["date"] = now
            commit = get_commit_id()
            print("commit is updated.")
            meta.update(commit)
    else:
        meta = make_meta()
    print("metadata", meta)

    if args.create_ref:
        ref_dic = {"result": None}
    else:
        if os.path.isfile(ref_file):
            with open(ref_file) as f:
                ref_dic = json.load(f)
        else:
            ref_dic = {"result": None}

    print(args, args.set)
    exe_list = exe_dic[args.set]
    exeutil = ExeUtil(exe_list)
    exeutil.show_exist()

    result_dic = OrderedDict()
    dump_result_dic = OrderedDict()
    chk_dic = OrderedDict()

    kkrtype = os.path.split(os.getcwd())[-1]
    akaikkr_path = os.path.join(args.program_path, kkrtype, akaikkr_exe)
    fmg_path = os.path.join(args.program_path, kkrtype, fmg_exe)
    print("executable files")
    print(akaikkr_path)
    print(fmg_path)
    print()
    if not os.path.isfile(akaikkr_path) or not os.path.isfile(fmg_path):
        print(f"failed to find {akaikkr_path} or {akaikkr_path}")
        sys.exit(100)

    prog = {"specx": akaikkr_path, "fmg": fmg_path, "args": args}
    execute_postscript = not args.no_postscript

    for exe_ in exe_list:
        key, result = exe_(prog, displc=displc,
                           execute_postscript=execute_postscript)
        result_dic[key] = result
        dump_result_dic[key] = deepcopy(result)
        chk = check_go(key, result_dic, ref_dic["result"])
        chk_dic[key] = chk

    print()
    print("================")
    print()
    if args.create_ref:
        with open(ref_file, "w") as f:
            json.dump({"meta": meta, "result": dump_result_dic}, f, indent=1)
        print("saved to", ref_file)
    else:

        chk_array = chk_dic_success_all(chk_dic)
        if np.all(chk_array):
            print("ALL TESTS PASSED.")
            print()
        else:
            print("TEST FAILED.")
            for chk_key, chk_value in chk_dic.items():
                if not chk_value:
                    print(chk_key)
            print()

        print("SHORT SUMMARY")
        resultutil = ResultUtil(chk_dic)
        resultutil.show()
        show_chk_legend()

        print("# reference:")
        for key, value in ref_dic["meta"].items():
            print("{}: {}".format(key, value))

        result_file = "result.json"
        with open(result_file, "w") as f:
            json.dump({"meta": meta, "result": dump_result_dic}, f, indent=1)
            print()
            print("# result saved to", result_file)
