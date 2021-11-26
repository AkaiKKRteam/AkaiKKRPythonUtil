# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import os
import subprocess
from pymatgen.core import Element, Structure
from pyakaikkr import AkaikkrJob
from pyakaikkr import DsplcMaker
from pyakaikkr.CompareCifKkr import CompareCifKkr
import pandas as pd
import numpy as np


def get_kkr_struc_from_cif(ciffilepath: str, akaikkr_exe: str, displc: bool,
                           use_bravais=True, remove_temperaryfiles=True,
                           directory="temporary") -> dict:
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
        parent_directory (str, optional): working directory.

    Returns:
        dict: kkr structure parameters on success.
    """
    comp = CompareCifKkr(ciffilepath, akaikkr_exe, displc=displc,
                         parent_directory=directory)
    result = comp.convert_and_compare(use_bravais=use_bravais)
    struc_param = None
    if result == comp.SUCCESS:
        struc_param = comp.get_structure_param(
            remove_temperaryfiles=remove_temperaryfiles)
        kkr_struc = comp.kkr_struc
    else:
        print("failed to convert the cif file")
        print("msg=", comp.msg)
        print("result=", result)
    try:
        os.rmdir(comp.parent_directory)
    except OSError:
        # ignore errors on removing output directory
        pass
    return struc_param, kkr_struc


def param_add_dsplc(param: dict, kkr_struc: Structure,
                    debye_temp: float, temp: float,
                    mxl: int = 3):
    """add dsplc to param

    Args:
        param (dict): kkr parameter
        kkr_struc (Structure): kkr structure
        debye_temp (float): Debye temperature
        temp (float): temperature
        mxl (int, optional): maximum value of l. Defaults to 3.

    Raises:
        ValueError: if it contans alloy.

    Returns:
        dict: kkr parameter
    """

    print("debug, param", param)
    brvtyp = param["brvtyp"]
    dsplcmaker = DsplcMaker(kkr_struc, param)

    dsplc_list = []
    ncmp_list = []
    mxl_list = []
    anclr_list = []
    conc_list = []
    for _Zlist, _conclist, _mxl in zip(param["anclr"], param["conc"], param["mxl"]):
        if len(_Zlist) == 1:
            pass
        else:
            raise ValueError("len(param[anclr]) must be 1")
        z = _Zlist[0]
        conc = _conclist[0]
        _mxl = mxl  # discard _mxl and use a new value
        elm = Element("H").from_Z(z)
        m = elm.atomic_mass
        displace = dsplcmaker.get_avr_debye(
            m, debye_temp, temp)
        if brvtyp == "bcc" or brvtyp == "fcc":
            dsplc = dsplcmaker.get_displc_as_cartesian(
                cart_len=displace,
                xyz=["x", "0", "0"],
                alen=param["a"],
            )
            ndsplc = len(dsplc)
        elif brvtyp == "hcp":
            dsplc = []
            if False:
                # if exact operations are known.
                d = displace/param["a"]
                dsplc1 = []
                dsplc1.append([d, 0, 0])
                dsplc1.append([-d, 0, 0])
                angle = np.pi/3.0
                dsplc1.append([d*np.cos(angle), d*np.sin(angle), 0])
                dsplc1.append([-d*np.cos(angle), -d*np.sin(angle), 0])
                angle = np.pi/3.0*2.0
                dsplc1.append([d*np.cos(angle),
                              d*np.sin(angle), 0])
                dsplc1.append(
                    [-d*np.cos(np.pi/3.0*2.0), -d*np.sin(np.pi/3*2.0), 0])
                dsplc1.append([0, 0, d])
                dsplc1.append([0, 0, -d])
                dsplc.extend(dsplc1)
            else:
                dsplc = dsplcmaker.get_displc_as_cartesian(
                    cart_len=displace, alen=param["a"],
                    xyz=[["x", "0", "0"],["0", "0", "z"]],  frac=True)

            print("dsplc", dsplc)

            ndsplc = len(dsplc)
        else:
            raise ValueError("unknown brvtyp={}".format(brvtyp))

        dsplc_list.append(dsplc)
        ncmp_list.append(ndsplc)
        mxl_list.append(_mxl)
        anclr_list.append([z]*ndsplc)
        conc_list.append([conc]*ndsplc)

    param["displc"] = dsplc_list
    param["ncmp"] = ncmp_list
    param["mxl"] = mxl_list
    param["anclr"] = anclr_list
    param["conc"] = conc_list
    print("  temperature: ", temp, "(K)")
    print("  mean square displaement: ", displace, "(bohr)")
    print("debug param=", param)
    return param


def load_debye_prop(filename="debye.csv"):
    df = pd.read_csv("debye.csv").set_index("elm")
    _df = df[["structure", "a(Ang.)", "c(Ang.)",
              "DebyeT(0K)", "ciffile", "magtyp"]]
    _df.dropna(how="any", axis=0, inplace=True)
    d = _df.to_dict()
    structure_dic = d["structure"]
    debye_temp_dic = d["DebyeT(0K)"]
    lattice_a_dic = d["a(Ang.)"]
    lattice_c_dic = d["c(Ang.)"]
    ciffile_dic = d["ciffile"]
    magtyp_dic = d["magtyp"]

    # change format
    lattice_const_dic = {}
    for elm in lattice_a_dic.keys():
        lattice_const_dic[elm] = [lattice_a_dic[elm], lattice_c_dic[elm]]
    del lattice_a_dic, lattice_c_dic

    return structure_dic, debye_temp_dic, lattice_const_dic, ciffile_dic, magtyp_dic


def main(akaikkr_exe):
    # target list
    target_list = ["Co"]
    temp_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

    structure_dic, debye_temp_dic, lattice_const_dic, ciffile_dic, magtyp_dic = load_debye_prop()

    # target loop
    for target in target_list:

        displc = True
        ciffilepath = ciffile_dic[target]
        struc_param, kkr_struc = get_kkr_struc_from_cif(
            ciffilepath, akaikkr_exe, displc)
        print(kkr_struc)

        elm = Element(target)
        z = elm.Z
        m = float(elm.atomic_mass)
        print("#", target)
        print(" ", structure_dic[target], lattice_const_dic[target], "(ang.)")

        debye_temp = debye_temp_dic[target]
        print("  atomic mass: ", m, " debye temp: ", debye_temp)

        # temperature loop
        for nrun, temp in enumerate(temp_list):
            directory = target+"/"+str(temp)
            os.makedirs(directory, exist_ok=True)
            job = AkaikkrJob(directory)
            param = job.default
            param["magtyp"] = magtyp_dic[target]

            param["sdftyp"] = "mjwasa"
            param["record"] = "init"
            param["edelt"] = 1e-3
            param["ewidth"] = 1.0
            param["bzqlty"] = 10
            param["pmix"] = 0.02

            param.update(struc_param)
            param["a"] = lattice_const_dic[target][0]/0.529177
            param["c/a"] = lattice_const_dic[target][1] / \
                lattice_const_dic[target][0]

            param = param_add_dsplc(param, kkr_struc, debye_temp, temp)

            option_param = {"cpaitr_show": True}
            param["option"] = option_param

            print("  directory: ", directory)

            inputcard_go = "inputcard_go"
            output_go = "out_go.log"
            if nrun == 0:
                # rough scf cal. for 1st concentrtion
                print("rough scf")

                job.make_inputcard(param,  inputcard_go)
                try:
                    job.run(akaikkr_exe, inputcard_go, output_go)
                except KKRFailedExecutionError as err:
                    print(err)
                    sys.exit(10)
            else:
                # simply copy converged potential file (scf tight)
                if True:
                    job.copy_potential(param["potentialfile"], param["potentialfile"],
                                       prev_dir, directory)
                else:
                    cmd = "cp {} {}".format(
                        prev_dir+"/"+param["potentialfile"], directory+"/"+param["potentialfile"])
                    subprocess.call(cmd, shell=True)

            # tight scf cal.
            param["record"] = "2nd"
            param["bzqlty"] = 20
            param["edelt"] = 1e-4
            print("tight scf")
            job.make_inputcard(param, inputcard_go)
            try:
                job.run(akaikkr_exe, inputcard_go, output_go)
            except KKRFailedExecutionError as err:
                print(err)
                sys.exit(10)

            if not job.check_convergence_go(output_go):
                print("failed to get converged, but continue.")
                continue

            # dos cal.
            print("dos")
            param["go"] = "dos"
            param["record"] = "2nd"
            param["bzqlty"] = 20
            param["edelt"] = 1e-4
            job.make_inputcard(param, "inputcard_dos")
            try:
                job.run(akaikkr_exe, "inputcard_dos", "out_dos.log")
            except KKRFailedExecutionError as err:
                print(err)
                sys.exit(10)

            # dos cal.
            print("dos")
            param["go"] = " cnd"
            param["record"] = "2nd"
            param["bzqlty"] = 20
            param["edelt"] = 1e-8
            param["ewidth"] = 0.01
            job.make_inputcard(param,  "inputcard_cnd")
            try:
                job.run(akaikkr_exe, "inputcard_cnd", "out_cnd.log")
            except KKRFailedExecutionError as err:
                print(err)
                sys.exit(10)

            prev_dir = directory


main(akaikkr_exe="/home/kino/kino/kit/AkaiKKRprogram.current.gfortran/akaikkr_cnd/specx")
