# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

#!/usr/bin/env python
# coding: utf-8


import json

import hashlib
import os
from collections import OrderedDict
import numpy as np


_PREFIX_ = os.getcwd()


def get_ciffiles(datacatalog="datacatalog/materiallibrary_smallset.json",
                 prefix=_PREFIX_):

    if True:
        with open(os.path.join(prefix, datacatalog)) as f:
            data = json.load(f)
    else:
        with open("datacatalog/materiallibrary_content.json") as f:
            data = json.load(f)
    return data


def cif2poscar(filename, primitive=True, directory=None,
               outputfilename=None, dry=False, fmt="poscar"):
    # input
    cirparser = CifParser(filename)
    struc = cirparser.get_structures(primitive=primitive)[0]

    # output
    if directory is None:
        hashstr = hashlib.md5(filename.encode()).hexdigest()
        directory = hashstr

    if not dry:
        os.makedirs(directory, exist_ok=True)
    if outputfilename is None:
        if fmt == "poscar":
            outputfilename = "poscar.vasp"
        elif fmt == "cif":
            outputfilename = "struc.cif"

    filepath = os.path.join(directory, outputfilename)
    if not dry:
        struc.to(fmt, filepath,)
    print("save to", filepath)
    return filepath


def read_poscar(filename, primitive=True):
    # input
    struc = Structure.from_file(filename)

    # analysis as it is
    spa = SpacegroupAnalyzer(struc)
    spainfo = (spa.get_space_group_symbol(), spa.get_space_group_number())

    # analysis as standard structure
    if primitive:
        stand_struc = spa.get_primitive_standard_structure()
    else:
        stand_struc = spa.get_conventional_standard_structure()
    spg = SpacegroupAnalyzer(stand_struc)
    spginfo = (spg.get_space_group_symbol(), spg.get_space_group_number())
    print("spa", spainfo, len(struc.sites),
          "spg", spginfo, len(stand_struc.sites))


def sep_filename(filename):
    s = os.path.split(filename)
    directory = s[0]
    filename = s[-1]
    t = os.path.splitext(filename)
    filename_base = t[0]
    filename_ext = t[-1]
    return {"directory": directory, "filename": filename,
            "base": filename_base, "ext": filename_ext}


def make_p1structures(datacatalog,
                      datapath_prefix,
                      newdirectory,
                      outputjson,
                      fmt="cif"):

    files = get_ciffiles(datacatalog=datacatalog)

    for file in files:
        filename = os.path.join(datapath_prefix,file["filename"])
        print(filename)
        filename_info = sep_filename(filename)
        if fmt == "cif":
            outputfilename = filename_info["base"]+".cif"
        elif fmt=="vasp" or fmt=="poscar":
            outputfilename = filename_info["base"]+".vasp"
        else:
            raise ValueError("unknown fmt={}".format(fmt))
        try:
            cif2poscar(filename, primitive=True, directory=newdirectory,
                       outputfilename=outputfilename,
                       dry=False)
        except ValueError as err:
            continue
        newfilename = os.path.join(newdirectory, outputfilename)
        if fmt == "cif":
            file["filename-p1-cif"] = newfilename
        elif fmt=="vasp" or fmt=="poscar":
            file["filename-poscar"] = newfilename

    with open(outputjson, "w") as f:
        json.dump(files, f, indent=1)

    return file


if __name__ == "__main__":
    def main(datapath_prefix):
        originaldatacatalog = "datacatalog/materiallibrary_smallset.json"
        newdirectory = "MaterialsLibrary-P1-smallset"
        outputdatacatalog = "datacatalog/materiallibrary_smallset_p1.json"
        if True:
            files = make_p1structures(datacatalog=originaldatacatalog,
                                    datapath_prefix=datapath_prefix,
                                    newdirectory=newdirectory,
                                    outputjson=outputdatacatalog,
                                    fmt="poscar")

        if False:
            with open(os.path.join(outputdatacatalog)) as f:
                files = json.load(f)

            # confirm that they are be read.
            for file in files:
                if "filename-poscar" in file:
                    read_poscar(file["filename-poscar"],)

    homedir = os.path.expanduser("~")
    datapath_prefix = os.path.join(homedir,"kino/kit/MaterialsLibrary")
    main(datapath_prefix)