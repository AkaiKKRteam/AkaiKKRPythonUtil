# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from pymatgen.core.periodic_table import Element
from collections import OrderedDict
import matplotlib.pyplot as plt
import re
import numpy as np


if False:
    # Awk can be read, but an error occurs in plotting

    def _load_spc_format_type_c(lines_iter):
        line = next(lines_iter)
        s = line.strip().split()
        emin = float(s[1])
        emax = float(s[2])
        ne = int(s[3])
        nhighsymkp = int(s[4])

        line = next(lines_iter)
        line = line.strip()
        kstr_list = re.split("  +", line)
        kcrt = OrderedDict()
        for kstr in kstr_list[1:]:
            print("kstr", kstr)
            s = kstr.split()
            idx = int(s[0])
            name = " ".join(s[1:]).replace("'", "")
            kcrt[idx] = name

        # omit checcking nhighsymkp == len(kcrt.keys())

        while True:
            line = next(lines_iter)
            if line.startswith("### end of header"):
                break

        Awk = []
        while True:
            try:
                line = next(lines_iter)
            except StopIteration:
                break
            line = line.strip()
            s = re.split(" +", line)
            Aw = list(map(float, s))
            Awk.append(Aw)

        Awk = np.array(Awk)
        return emin, emax, ne, nhighsymkp,  kcrt, Awk
        # reutrn kcrt, Awk, kdist, energy, kpath

    def _load_spc(filepath, format_type):
        with open(filepath) as f:
            lines = f.read().splitlines()

            lines_iter = iter(lines)
            line = next(lines_iter)
            s = line.strip().split()
            internal_format_type = s[-1]

            if format_type == "(c)":
                if internal_format_type == format_type:
                    emin, emax, ne, kcrt, highsymkp, Awk = _load_spc_format_type_c(
                        lines_iter)
            else:
                print("format type {} is not supported.", format_type)
                raise ValueError

        return emin, emax, ne, kcrt, Awk


def _load_spc_format_type_a_Awk(lines_iter: iter):
    """load A(w,k) in the spc format type a from iterator

    Args:
        lines_iter (iter): iterator of kkr outputs

    Returns:
        np.ndarray, np.ndarray, np.ndarray: Awk, kdist, energy
    """
    Awk = []
    Aw = []
    for x in lines_iter:
        if x == "" and len(Aw) > 0:
            Awk.append(Aw)
            Aw = []
        if x != "":
            y = re.split(" +", x.strip())
            y = list(map(float, y))
            Aw.append(y)
    if len(Aw) > 0:
        Awk.append(Aw)

    kdist = []
    for Aw in Awk:
        kdist.append(Aw[0][0])
    energy = []
    for Aw in Awk[0]:
        energy.append(Aw[1])
    Awk2 = []
    for Aw in Awk:
        wlist = []
        for w in Aw:
            wlist.append(w[2])
        Awk2.append(wlist)
    return np.array(Awk2), np.array(kdist), np.array(energy)


def _load_spc_format_type_a(filepath: str):
    """load A(w,k) in the spc format type a

    Args:
        filepath (str): output filename

    Returns:
        np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray: kcrt, Awk, kdist, energy, kpath
    """
    with open(filepath) as f:
        lines = f.readlines()
    lines2 = []
    for line in lines:
        lines2.append(line.strip())
    lines = lines2

    lines_iter = iter(lines)
    line = next(lines_iter)
    s = line.strip().split()
    internal_format_type = s[-1]

    if internal_format_type == "(a)":
        line = next(lines_iter)
        if False:
            s = line.strip().split()
            emin = float(s[1])
            emax = float(s[2])
            ne = int(s[3])
            nhighsymkp = int(s[4])

        line = next(lines_iter)
        line = line.strip()
        kstr_list = re.split("   +", line)
        kpath = OrderedDict()
        kcrt = []
        for kstr in kstr_list[1:]:
            s = kstr.split()
            idx = int(s[0])
            idx -= 1  # index convert from fortran to Python
            name = " ".join(s[1:]).replace("'", "")
            kpath[idx] = name
            kcrt.append(idx)
        kcrt = np.array(kcrt)
    elif internal_format_type == "(a1)":  # short format without kpoints
        line = next(lines_iter)
        line = line.strip()
        kstr_list = line.split()
        kpath = OrderedDict()
        kcrt = []
        kpoints = None
        for kstr in kstr_list[2:]:
            s = kstr.split()
            idx = int(s[0])
            idx -= 1  # index convert from fortran to Python
            name = None
            kpath[idx] = name
            kcrt.append(idx)
        kcrt = np.array(kcrt)

    while True:
        line = next(lines_iter)
        if line.startswith("### end of header"):
            break
    line = next(lines_iter)
    Awk, kdist, energy = _load_spc_format_type_a_Awk(lines_iter)

    return kcrt, Awk, kdist, energy, kpath


class AwkReader:
    """read A(w,k)
    """

    def __init__(self, filepath, fmt=1):
        """initialization routine

        Args:
            filepath (str): filepath of Akaikkr output filename run as go="spc*"
        """

        if fmt == 1:
            self.kcrt, self.Awk, self.kdist, self.energy, self.kpath = _load_spc_format_type_a(
                filepath)
        else:
            print("unsuported format {}".format(fmt))
            raise ValueError
