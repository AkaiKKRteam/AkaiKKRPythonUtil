# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import sys
import os
import json

if "test_script" not in sys.modules:
    from pyakaikkr.CompareCifKkr import CompareCifKkr


def get_kkr_struc_from_cif(ciffilepath: str, specx: str, displc: bool,
                           use_bravais=True, remove_temperaryfiles=True,
                           fmt="cif"):
    """get kkr structure parameter from cif file.

    If specx is akaikkr, then displc is set to False.
    If specx is akaikkr_cnd, then displc is set to True.

    If use_bravais is True, bcc,bcc,.., and a,c,b,alpha,beta,gaam are used. 
    If use_bravias is False, aux and lattice vectors iare used.

    rmt, lmax, ..., displc are set to be default values. 

    Args:
        ciffilepath (str): cif file path
        specx (str, optional): specx path. Defaults to str.
        displc (bool, optional): displc parameter is added. 
        use_bravais (bool, optional): use bravias lattice. Defaults to True.
        remove_temperaryfiles (bool, optional): delete temporary files on exit. Defaults to True.
        fmt (str, optional): file format. Defaults to "cif".
    Returns:
        dict: kkr structure parameters on success.
    """
    comp = CompareCifKkr(ciffilepath, specx, displc=displc, fmt=fmt)
    result = comp.convert_and_compare(use_bravais=use_bravais)
    struc_param = None
    if result == comp.SUCCESS:
        struc_param = comp.get_structure_param(
            remove_temperaryfiles=remove_temperaryfiles)
    else:
        print("failed to convert the cif file")
        print("msg=", comp.msg)
        print("result=", result)
        sys.exit(10)

    try:
        os.rmdir(comp.parent_directory)
    except OSError:
        # ignore errors on removing output directory
        pass
    return struc_param


if __name__ == "__main__":
    def main(path_prefix, ciffile_path, akaikkr_type="akaikkr",):
        """load data from path_prefix and convert to the PyAkaiKKR dict format.

        The output is printed to stdout.

        akaikkr_type can be akaikkr or akaikkr_cnd.

        Args:
            path_prefix (str): path prefix to AkaiKKR
            ciffile_path (st): the cif file name.
            akaikkr_type (str, optional): type of AkaiKKR. Defaults to "akaikkr".

        Raises:
            ValueError: unknown fmt.
        """
        if akaikkr_type == "akaikkr":
            displc = False
        elif akaikkr_type == "akaikkr_cnd":
            displc = True
        else:
            raise ValueError("unknown akaikkr_type={}".format(akaikkr_type))
        specx = os.path.join(path_prefix, akaikkr_type, "specx")

        use_bravais = True

        struc_param = get_kkr_struc_from_cif(ciffile_path, specx,
                                             use_bravais=use_bravais,
                                             displc=displc, fmt="cif")
        print()
        print("sturc_param")
        if False:
            for key, value in struc_param.items():
                print(key, value)
        else:
            print(json.dumps(struc_param))

    def define_and_get_parse():
        import argparse
        parser = argparse.ArgumentParser()

        parser.add_argument("--akaikkr", default= "kino/kit/AkaiKKRprogram.current.gfortran")
        parser.add_argument("--prefix", default= "kino/kit/MaterialsLibrary")
        parser.add_argument("--akaikkr_type", choices = ["akaikkr", "akaikkr_and"], default= "akaikkr")

        parser.add_argument("structure_file") 
        # e.g. ="kino/kit/MaterialsLibrary/MaterialsLibrary/AtomWorkData/small_sites/made_by_kino/Co_P63mmc.cif"

        args = parser.parse_args()
        return args

    homedir = os.path.expanduser("~")
    args = define_and_get_parse()
    main(os.path.join(homedir, args.akaikkr),
         os.path.join(homedir, args.prefix, args.structure_file), 
         args.akaikkr_type)
