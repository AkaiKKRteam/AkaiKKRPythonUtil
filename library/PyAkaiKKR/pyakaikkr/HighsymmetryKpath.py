# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from collections import OrderedDict
from pymatgen.core import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pymatgen.symmetry import kpath
import json
import numpy as np


def _fix_kname(name):
    """fix kname to the plottable one
    only for Gamma

    Args:
        name (str): kpath name

    Returns:
        str: fixed kpath name
    """    
    if name=="GAMMA":
        name = "\\Gamma"
    return name

class HighSymKPath:
    def __init__(self, structure, path_type=None, klabel_filename="kpath.json"):
        """initilization routine

        path_type is path_type of HighSymmKpath

        Args:
            structure (Structure): pymatgen.core.Structure
            path_type (str): path_type. Defaults to None.
            klabel_filename ((str,io.TextIOBase), optional): filename with klabel. Defaults to "kpath.json".
        """
        if path_type is None:
            print("kpath is made by HighSymmKpath")
            _kpath = HighSymmKpath(structure)
        elif path_type == "seekpath" or path_type == "hinuma":
            print("kpath is made by kpath.KPathSeek")
            _kpath = kpath.KPathSeek(structure)
        else:
            print("kpath is made by HighSymmKpath, path_type", path_type)
            _kpath = HighSymmKpath(structure, path_type=path_type)

        if False:
            # klabel will be saved in make_akaikkr_lines
            if klabel_filename is not None:
                with open(klabel_filename, "w") as f:
                    content = OrderedDict()
                    for key, value in _kpath.kpath["kpoints"].items():
                        value = np.array(value) # value is list or np.array
                        if key=="GAMMA":
                            key = "\\Gamma"
                        content[key] = value.tolist()
                    print("debug: content",content)
                    json.dump(content, f, indent=2)
                    print("klabel saved to", klabel_filename)
        self.kpath = _kpath
        self.klabel_filename = klabel_filename

    def make_akaikkr_lines(self, fmt=3, nk=250, nk_each=30, first_connected_kpath=True):
        """generate high symmetry kpath in the Akaikkr type(c) format.

        Args:
            fmt (int, optional): format type 3=type(c). Defaults to 3.
            nk (int, optional): the number of total kpoints for type(c). Defaults to 250.
            nk_each (int, optional): the number of total kpoints for type(b). Defaults to 30.
            first_connected_kpath (bool, optional): use only the first k path. Defaults to True.

        Raises:
            ValueError: uknown format type
        Returns:
            [str]: kkr input format
        """        

        _kpath = self.kpath
        klabel_filename = self.klabel_filename
        content = []
        lines = []

        if fmt == 2:
            for vpath in _kpath.kpath["path"]:
                s = ["c --- "]
                s.extend(vpath)
                lines.append(" ".join(s))
                acontent = []
                for ik, k in enumerate(vpath):
                    s = ["c --- ", k]
                    lines.append(" ".join(s))
                    v = _kpath.kpath["kpoints"][k]
                    v = np.array(v)
                    s = list(map(str, v))
                    if ik == 0:
                        s.append(str(0))
                    else:
                        s.append(str(nk_each))
                    lines.append(" ".join(s))
                    k = _fix_kname(k)
                    acontent.append({k: v.tolist()})
                content.append(acontent)
                if first_connected_kpath:
                    break
        elif fmt == 3:
            lines = [str(nk)]
            # use only the first path by "break" inside for
            for vpath in _kpath.kpath["path"]:
                s = ["c --- "]
                s.extend(vpath)
                lines.append(" ".join(s))
                for k in vpath:
                    s = ["c --- ", k]
                    lines.append(" ".join(s))
                    v = _kpath.kpath["kpoints"][k]
                    v = np.array(v)
                    s = list(map(str, v))
                    lines.append(" ".join(s))
                    k = _fix_kname(k)
                    content.append({k: v.tolist()})
                content = [content]
                break
        else:
            raise ValueError("unknown fmt={}".format(fmt))

        
        if klabel_filename is not None:
            import io
            if isinstance(klabel_filename, str):
                with open(klabel_filename, "w") as f:
                    json.dump({"kpath": content}, f)
                # print("  klabel saved to", klabel_filename)
            elif isinstance(klabel_filename, io.TextIOBase):
                json.dump({"kpath": content}, klabel_filename)

        if False:
            klist = []
            for v in content:
                k = str(list(v.keys())[0])
                klist.append(k)

        return lines
