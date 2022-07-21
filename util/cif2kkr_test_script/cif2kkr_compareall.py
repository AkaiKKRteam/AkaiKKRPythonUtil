# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


import json
import sys
import os


"""
requires 
spglib                             >=1.16.2
pymatgen                           >=2022.0.14
"""

if "test_script" not in sys.modules:
    from pyakaikkr import AkaikkrJob
    from pyakaikkr.Error import *
    from pyakaikkr import StructureSpeciesConverter
    from pyakaikkr.CompareCifKkr import *


def get_ciffiles(prefix, dataset="smallset"):
    """load cif filenames.

    dataset name can be "smallset", "all", "smallset_p1", "atomwork".

    Args:
        prefix (str): directory name
        dataset (str, optional): data set name. Defaults to "smallset".

    Returns:
        output contains

        - str: file names
        - str: file format, cif or poscar for P1.
    """
    fmt = "cif"
    # from list

    # A small size cif files covering all the space group numbers
    if dataset == "smallset":
        with open(os.path.join(prefix, "datacatalog/materiallibrary_smallset.json")) as f:
            data = json.load(f)
        files = []
        for datum in data:
            files.append(os.path.join(prefix, datum["filename"]))
        print("datasize=", len(files))
    # atomwork and additional data
    elif dataset == "atomwork":
        with open(os.path.join(prefix, "datacatalog/atomworkdata_content.json")) as f:
            data = json.load(f)
        files = []
        for datum in data:
            files.append(os.path.join(prefix, datum["filename"]))
        print("datasize=", len(files))
 
    # all the materiallibrary
    elif dataset == "all":
        with open(os.path.join(prefix, "datacatalog/materiallibrary_content.json")) as f:
            data = json.load(f)
        files = []
        for datum in data:
            files.append(os.path.join(prefix, datum["filename"]))
        print("datasize=", len(files))
    # P1 cif generatef rom  materiallibrary_smallset
    elif dataset == "smallset_p1":
        datacatalogfile = os.path.join(
            prefix, "datacatalog/materialsLibrary_smallset_p1.json")
        with open(datacatalogfile) as f:
            data = json.load(f)

        files = []
        key = "filename-poscar"
        for datum in data:
            if key in datum:
                files.append(os.path.join(prefix, datum[key]))
        fmt = "poscar"
        print("datasize=", len(files))

    elif dataset == "set1":
        # spacial files which failed at some stage.
        # set 1
        files = [
            "MaterialsLibrary/chalcogenides3/Rb5AgTi6Se27 (153K,P31c).cif",
            "MaterialsLibrary/chalcogenides1/(NH4)2{Ni(en)2}3{Re6Te8.cif",
        ]
    elif dataset == "set2":
        # set 2
        files = ["MaterialsLibrary/chalcogenides2/LaTe2 (Super,Pc).cif",
                 "MaterialsLibrary/chalcogenides3/RbTaCu2Te4 (P21cn).cif",
                 "MaterialsLibrary/chalcogenides2/H2S (HP,LT,IV',3.8GPa,27K.cif"]
    elif dataset == "set3":
        files = [
            "MaterialsLibrary/chalcogenides2/Ho8Se14.7 (295K,Amm2).cif"]
    elif dataset == "set4":
        files = [
            "MaterialsLibrary/chalcogenides3/Se3.02S4.98 (P2&c).cif"]
        # have Q
    elif dataset == "set5":
        files = [
            "MaterialsLibrary/chalcogenides2/In11Sn5.5S22 (P2&m).cif"]
    elif dataset == "set6":
        files = [
            "MaterialsLibrary/chalcogenides2/K2Se5 (183K,P212121).cif"]
    elif dataset == "set7":
        files = [
            "MaterialsLibrary/AtomWorkData/spcno_series/spcno_023_BaAG1.5Cu0.5SnS4.cif"]
    elif dataset == "set8":
        files = [
            'MaterialsLibrary/chalcogenides3/ReS2,Rheniite (II,P-1).cif']
    elif dataset == "set9":
        files = ["MaterialsLibrary/AtomWorkData/spcno_series/spcno_172_NaCo.cif"]
    elif dataset == "set10":
        files = ["MaterialsLibrary/Intermetallics/Co0.60Ni0.40Sn2 (Aba2).cif"]
    elif dataset == "set11":
        files = [
            "MaterialsLibrary/AtomWorkData/spcno_series/spcno_177_Al4As12Ca3H30MgNa4O57.cif"]
    elif dataset == "set12":
        files = ["MaterialsLibrary/Intermetallics/Mn5Si2 (P41212).cif"]
        # space group is different.
    elif dataset == "set13":
        files = ["MaterialsLibrary/Appendix/Ba1.8Ca2.2Si6N10O (P213).cif"]
    elif dataset == "set14":
        files = ["MaterialsLibrary/Appendix/(Hg0.82Cu0.18)Ba2Ca2Cu3O8.cif"]
    elif dataset == "set15":
        files = ["MaterialsLibrary/Appendix/AsTe0.5O2 (P3).cif"]
    elif dataset == "set16":
        files = ["MaterialsLibrary/Appendix/HBr (Ib,116K,Pa3).cif"]
    elif dataset == "set17":
        files = ["MaterialsLibrary/Intermetallics/CeMg1.03 (P63&mmc).cif"]
    elif dataset == "set18":
        files = ["MaterialsLibrary/Intermetallics/Mn5Si2 (P41212).cif"]
    elif dataset == "set19":
        files = ["MaterialsLibrary/Intermetallics/Np (Beta,P4212).cif"]
    elif dataset == "set20":
        files = ["MaterialsLibrary/Intermetallics/Zr2Ni2In (P42&mnm).cif"]
    elif dataset == "set21":
        files = ["MaterialsLibrary/Pnictides/GdCuAs1.15P0.85 (Pmmn).cif"]
    elif dataset == "set22":
        files = ["MaterialsLibrary/Pnictides/Li6.45Mn3As4 (Cmma).cif"]
    elif dataset == "set23":
        files = ["MaterialsLibrary/chalcogenides1/(BPh4)(PEt3)5Re6Se8(CN).cif"]
    elif dataset == "set24":
        files = [
            "MaterialsLibrary/chalcogenides1/[Cu3Mo(O)S3(C5H4NS)(PPh3)3.cif"]
    elif dataset == "set25":
        files = [
            "MaterialsLibrary/AtomWorkData/spcno_series/spcno_042_Bi5FeO15Ti3.cif"]
    elif dataset == "set26":
        files = ["MaterialsLibrary-P1-smallset/Tb16Ni36P22 (P-6m2).cif",
                 "MaterialsLibrary-P1-smallset/Li5Ag0.094Hf3S8 (P4332).cif", ]
    elif dataset == "set27":
        files = ["MaterialsLibrary-P1-smallset/Tb16Ni36P22 (P-6m2).vasp",
                 "MaterialsLibrary-P1-smallset/Li5Ag0.094Hf3S8 (P4332).vasp",
                 ]
        fmt = "poscar"
    return files, fmt


def append_result_to_json(results, resultfile="result.json"):
    """append result to the json file.

    Args:
        results (dict): result of comparison
        resultfile (str, optional): output filename. Defaults to "result.json".
    """
    with open(resultfile, "w") as f:
        json.dump(results, f, indent=1)


def load_results(resultfile="result.json"):
    """load dict from the resultfile. 

    Args:
        resultfile (str, optional): result filename. Defaults to "result.json".

    Returns:
        dict: result dict
    """
    results = []
    if os.path.isfile(resultfile):
        with open(resultfile, "r") as f:
            results = json.load(f)
    return results


if __name__ == "__main__":
    def main(path_prefix, datapath_prefix, akaikkr_type="akaikkr", dataset="all"):
        """main routine of this program.

        It compares the structures from the input file and AkaiKKR output file.

        akaikkr_type can be "akaikkr" or "akaikkr_cnd" because their input formats are different.

        Args:
            path_prefix (str): directory path for akaikkr specx program.
            datapath_prefix (str): data path for structures
            akaikkr_type (str, optional): akaiKKR program type. Defaults to "akaikkr".

        Raises:
            ValueError: _description_
        """
        if akaikkr_type == "akaikkr":
            displc = False
        elif akaikkr_type == "akaikkr_cnd":
            displc = True
        else:
            raise ValueError("unknown akaikkr_type={}".format(akaikkr_type))

        specx = os.path.join(path_prefix, "{}/specx".format(akaikkr_type))

        use_bravais = True

        ciffiles, fmt = get_ciffiles(datapath_prefix, dataset)

        special_ciffiles = []

        results = load_results()
        processed_ciffiles = []
        if len(results) > 0:
            processed_ciffiles = [x["ciffile"] for x in results]

        logfile = "log.txt"

        # for ciffile in ciffiles:
        for ciffile in ciffiles:

            if ciffile in processed_ciffiles:
                # print("skip processed file", ciffile)
                continue

            if ciffile not in special_ciffiles:
                print("========================", ciffile)
                result = None
                comp = CompareCifKkr(ciffile, specx, displc=displc, fmt=fmt)

                print("===", comp.directory)

                result = comp.convert_and_compare(use_bravais=use_bravais)
                print("msg=", comp.msg)
                print("result=", result)

                if comp.input_struc is None:
                    result = comp.NO_STRUCTURE

                results.append({"directory": comp.directory,
                                "ciffile": ciffile,
                                "result": result})

                append_result_to_json(results)

                if result != comp.SUCCESS:
                    comp.log_msg(comp.msg, logfile)

                print()
        print()
        print("all done")

    def define_and_get_parse():
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument("--akaikkr", default= "kino/kit/AkaiKKRprogram.current.gfortran")
        parser.add_argument("--prefix", default= "kino/kit/MaterialsLibrary")
        parser.add_argument("--datatype", choices = ["smallset", "atomwork", "all", "smallset_p1"], default= "all")
        parser.add_argument("--akaikkr_type", choices = ["akaikkr", "akaikkr_cnd"], default= "akaikkr")

        args = parser.parse_args()
        return args

    homedir = os.path.expanduser("~")
    args = define_and_get_parse()

    main(os.path.join(homedir, args.akaikkr),
         os.path.join(homedir, args.prefix),
         args.akaikkr_type, args.datatype)
