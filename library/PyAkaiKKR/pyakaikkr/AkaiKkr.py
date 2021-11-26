# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from pymatgen.core import Structure
import subprocess
import os
import pandas as pd
import numpy as np
import shutil
from scipy import integrate

from pymatgen.core.periodic_table import Element
from .AwkReader import AwkReader
from .Cif2Kkr import ak_cif2kkrparam
from .Error import KKRValueAquisitionError, KKRFailedExecutionError
from .Unit import Unit

# from .DsplcMaker import DsplcMaker

_au2ang = Unit().length_au2ang


def _fix_uniq_atom_names(uniq_atom_names):
    """fix uniq atom names

    convert  Composition1, Composition2, ... to Element1, Element2, ...

    Args:
        uniq_atom_names (str): a list of Composition

    Returns:
        [str]: a list of Elements
    """
    """ fix uniq_atom_names for pymatgen read poscar"""
    names = []
    Z = 0
    for i in range(len(uniq_atom_names)):
        Z += 1
        if Z == 2 or Z == 10 or Z == 18 or Z == 36 or Z == 54 or Z == 86:  # skip noble gas
            Z += 1
        elm = Element.from_Z(Z)
        names.append(str(elm))
    return names


def _noalloy_atom_names(atom_names):
    """convert atom names, which may be Compostition, to Element names

    Args:
        atom_names ([Str]): a list of atom names

    Returns:
        [Str]: a list of atom names
    """
    _uniq_atom_names = list(set(atom_names))
    _changed_atom_names = _fix_uniq_atom_names(_uniq_atom_names)
    # make conversion dict
    conversion = {}
    for src, tgt in zip(_uniq_atom_names, _changed_atom_names):
        conversion[src] = tgt
    output_atom_names = []
    for src in atom_names:
        output_atom_names.append(conversion[src])
    return output_atom_names


def _analyze_float_format(line: str):
    """find the position of space
    by ind the location of space and the location of non space.

    Args:
        line (str): a line

    Returns:
        [(int,int)]: a list of [start potition, end position]
    """

    import re
    # search_str = r"(-?[0-9]+\.?[0-9]+)"
    search_str = r"[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?"
    rex = re.compile(search_str)
    keywords = rex.findall(line)
    span_list = []
    istart = 0
    for key in keywords:
        m = rex.search(line[istart:])
        span_list.append([istart, istart+m.end()])
        istart += m.end()
    return span_list


if True:
    def _analyze_pdos_format(keyword, data):
        """analyze pdos format
        Sometimes **** is included, or digits are collapsed.
        The format is obtained by using the row with the highest split number from the first 10 rows.

        Args:
            keyword (str): keyword string
            data ([str]): lines

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [int, int]: a list of [start:end] columns
        """

        icmp = 1
        form_list = []
        still_continue = True
        while still_continue:
            search_keyword = " {}{:2d}".format(keyword, icmp)
            data_iter = iter(data)
            while still_continue:
                try:
                    line = next(data_iter).rstrip()
                except StopIteration:
                    still_continue = False
                    break

                if search_keyword == line:
                    lines = []
                    # ****が入っていたり，や桁が潰れていることがある．
                    # 最初の10行から最もsplit数が大きな行を利用してformatを得る．
                    for i in range(10):
                        line = next(data_iter).rstrip()
                        lines.append(line)
                    # search max. number of line items
                    n_lines = [len(v.split()) for v in lines]
                    i_max = np.argmax(n_lines)
                    form = _analyze_float_format(lines[i_max])
                    form_list.append(form)
                    icmp += 1
                    break
        if len(form_list) == 0:
            raise KKRValueAquisitionError("failed to get pdos format")
        return form_list
else:
    def _analyze_pdos_format(keyword_list, data):
        if not isinstance(keyword_list, list):
            keyword_list = [keyword_list]
        check = [False for x in keyword_list]
        # cut the first PDOS line
        for ikey, keyword in enumerate(keyword_list):
            for line in data:
                if keyword in line:
                    check[ikey] = True
                if not line:
                    check[ikey] = False
                if check[ikey] and keyword not in line:
                    break
            if all(check):
                break

        span_list = _analyze_float_format(line)
        return span_list


def _cut_a_pdos_block(data_iter, search_keyword, value_positions):
    """Isolate the PDOS block

    accept that no search_keyword is found. It occurs in reading the down component in the case of nonmag.

    Args:
        data_iter (Iter): iterator of lines
        search_keyword (str): a string to find PDOS
        value_postiions ([int,int]): columns to read PDOS values.

    Returns:
        [list],[[list]]: energy, pdos values
    """
    while True:
        try:
            line = next(data_iter)
            line = line.rstrip()
        except StopIteration:
            break
        if line.startswith(search_keyword):  # found keyword
            # 空白行まで読み込む
            lines = []
            while True:
                line = next(data_iter)
                line = line.strip()
                if len(line) == 0:
                    break
                lines.append(line)
            # value_positonsに合わせてfloatに分ける．
            e = []  # energy has a common values for all pdos
            pdos = []
            for line in lines:
                values = [line[istart:iend]
                          for istart, iend in value_positions]
                e.append(float(values[0]))
                v2 = []
                for s in values[1:]:
                    if s.startswith("**********"):
                        v2.append("1000.0")
                    else:
                        v2.append(s)
                try:
                    v = list(map(float, v2))
                except ValueError:
                    # I want to write outfile here...
                    print("\nerror in cut_pdos_all")
                    return None, None
                pdos.append(v)
            if len(e) == 0:
                return None, None
            return e, pdos

    return None, None


def _cut_only_nn(_jijdf):
    """cut only the first columns if jij DataFrame

    Args:
        _jijdf (pd.DataFrame): Jij DataFrame

    Returns:
        pd.DataFrame: Jij only the first n.n. neighbor
    """
    _jijdf.reset_index(inplace=True)
    type1 = _jijdf["type1"].values
    type2 = _jijdf["type2"].values
    typepair = []
    for t1, t2 in zip(type1, type2):
        typepair.append("-".join([t1, t2]))
    _pairdf = pd.DataFrame({"pair": typepair})
    jijdf = pd.concat([_jijdf, _pairdf], axis=1)
    typepair = list(set(typepair))
    typepair.sort()
    df_list = []
    for pair in typepair:
        _df = jijdf.query("pair=='{}'".format(pair)).sort_values(
            by="distance").reset_index(drop=True)
        _dist = _df.loc[0, "distance"]
        df = _df[_df["distance"] == _dist]
        df_list.append(df)
    dfsmall = pd.concat(df_list, axis=0)
    del dfsmall["index"]
    jijnn = [list(dfsmall.columns)]
    v = dfsmall.values.tolist()
    jijnn.extend(v)
    return jijnn


class AkaikkrJob:
    """Akaikkr job class to generate input, to run, and to analyze output
    """

    def __init__(self, path_dir):
        """initialization routine

        Args:
            path_dir (str): directory to run
        """
        # path of running directory
        self.path_dir = path_dir
    # default parameters for AkaiKKR
        self.default = {
            "go": "go", "potentialfile": "pot.dat",
            "brvtyp": "bcc", "a": 0, "c/a": 1.0, "b/a": 1.0,
            "alpha": 90, "beta": 90, "gamma": 90,
            "edelt": 1e-3, "ewidth": 1.0,
            "reltyp": "sra", "sdftyp": "mjw", "magtyp": "nmag", "record": "2nd",
            "outtyp": "update", "bzqlty": 6, "maxitr": 200, "pmix": "0.02",
            "ntyp": 1,
            "rmt": [1.0],
            "field": [0.0],
            "mxl": [2],
            "type": ["Fe"],
            "ncmp": [1],
            "anclr": [[26]],
            "conc": [[100]],
            "natm": 1,
            "atmicx": [
                ["0.00a", "0.00b", "0.00c", "Fe"],
            ],
        }

    def _read(self, outfile):
        """read outfile

        Args:
            outfile (str): output filename

        Returns:
            list: newline splitted contents of the file
        """
        filepath = os.path.join(self.path_dir, outfile)
        if True:
            data = None
            with open(filepath) as f:
                data = f.read().splitlines()
            return data
        else:
            try:
                with open(filepath) as f:
                    data = f.read().splitlines()
                return data
            except FileNotFoundError as e:
                print(e)
                print("FileNotFoundError: filename=", filepath, e)
                return []
            except Exception as e:
                print(e)
                print("Exception: filename=", filepath, e)
                return []

    def read_structure(self, structurefile, fmt="cif", use_bravais=True, use_primitive=True,
                       cif_primitive=True):
        """read structure and output it in the kkr input format

        if use_bravais==True, the output uses fcc, bcc and so on.
        otherwise, aux is used.

        Args:
            structurefile (str): structure filename
            fmt (str, optional): format of the structure filename. Defaults to "cif".
            use_bravais (bool, optional): use bravais lattice. Defaults to True.
            use_primitive (bool, optional): use primitive lattice. Defaults to True.
            cif_primitive (bool, optional): read structure as primitive lattice. Defaults to True.

        Returns:
            dict: kkr structure parameter
        """
        if fmt == "cif" or fmt == "vasp" or fmt == "poscar":
            self.param = ak_cif2kkrparam(structurefile, use_bravais=use_bravais,
                                         use_primitive=use_primitive,
                                         cif_primitive=cif_primitive, fmt=fmt)
        return self.param

    def make_inputcard(self, dic, inputcard):
        """make inputcard from dic

        Args:
            dic (dict): parameter of akaikkr input
            inputcard (str): input filename
        """

        def from_keylist(dic, keylist):
            """make dic.keys line and dic.values line separetely

            Args:
                dic (dict): parameters
                keylist ([str]): a list of kkr parameters

            Returns:
                [str]: keys line and values line
            """
            s = []
            result = []
            head = ["#---"]
            for key in keylist:
                head.append(key)
                s.append(str(dic[key]))
            result.append(" ".join(head))
            result.append(" ".join(s))
            return result

        def make_lattice_vectors(dic):
            """make r1,r2,r3 and a section from dic

            Args:
                dic (dict): parameters

            Returns:
                [str]: a list of kkr lattice parameters
            """
            result = []
            for vec in ["r1", "r2", "r3"]:
                s = dic[vec]
                result.append(" ".join([str(i) for i in s]))
            result.append(str(dic["a"]))
            return result

        def make_type_card(dic, longformat=False):
            """make type card section from dic

            if displc is in the dic, also add displc parameters

            Args:
                dic (dict): kkr parameters

            Returns:
                [str]: lines in kkr input format
            """
            result = []
            nn = 0
            if "displc" in dic:
                for i, j, k, l, m, _ in zip(dic["rmt"], dic["field"],
                                            dic["mxl"], dic["type"], dic["ncmp"],
                                            dic["displc"]):
                    if longformat:
                        result.append(
                            "#--- type ncmp rmt field mxl [displc anclr conc]")
                        xlist = [" ", l, str(m), str(i), str(j), str(k)]
                        for n in range(dic["ncmp"][nn]):
                            dispval = list(map(str, dic["displc"][nn][n]))
                            dispstr = " ".join(dispval)
                            xlist.extend(
                                [dispstr, str(dic["anclr"][nn][n]), str(dic["conc"][nn][n])])
                        result.append(" ".join(xlist))
                    else:
                        result.append("#--- type ncmp")
                        result.append(" ".join([" ", l, str(m)]))
                        result.append("#- rmt field mxl")
                        result.append(" ".join([str(i), str(j), str(k)]))
                        result.append("#- anclr conc")
                        for n in range(dic["ncmp"][nn]):
                            dispval = list(map(str, dic["displc"][nn][n]))
                            dispstr = " ".join(dispval)
                            result.append(
                                " ".join([dispstr, str(dic["anclr"][nn][n]), str(dic["conc"][nn][n])]))
                    nn += 1

            else:
                for i, j, k, l, m, in zip(dic["rmt"], dic["field"], dic["mxl"], dic["type"], dic["ncmp"]):
                    if longformat:
                        result.append(
                            "#--- type ncmp rmt field mxl [anclr conc]")
                        xlist = [" ", l, str(m), str(i), str(j), str(k)]
                        for n in range(dic["ncmp"][nn]):
                            xlist.extend(
                                [str(dic["anclr"][nn][n]), str(dic["conc"][nn][n])])
                        result.append(" ".join(xlist))
                    else:
                        result.append("#--- type ncmp")
                        result.append(" ".join([" ", l, str(m)]))
                        result.append("#- rmt field mxl")
                        result.append(" ".join([str(i), str(j), str(k)]))
                        result.append("#- anclr conc")
                        for n in range(dic["ncmp"][nn]):
                            # result.append(" ".join([" ".join(map(str,d[n])),str(dic["anclr"][n]),str(dic["conc"][n])]))
                            # result.append(" ".join(["0.0 0.0 0.0",str(dic["anclr"][n]),str(dic["conc"][n])]))
                            result.append(
                                " ".join([str(dic["anclr"][nn][n]), str(dic["conc"][nn][n])]))
                    nn += 1
            return result

        def make_atom_card(dic):
            """make atom card section from dic

            Args:
                dic (dict): kkr parameters

            Returns:
                [str]: lines in the kkr input format
            """
            result = []
            result.append("#--- atmicx atmtyp")
            for n in range(dic["natm"]):
                # s=[" ".join(dic["atmicx"][n]),dic["atmtyp"][n]]
                s = dic["atmicx"][n]
                result.append(" ".join(s))
            return result

        def make_option_card(dic):
            """make option section from dict.
                The option section starts with begin_option and ends with end_option
            Args:
                dic (dict): kkr parameters

            Returns:
                [str]: lines in the kkr input format
            """
            result = []
            if "option" in dic:
                options = dic["option"]
                result.append("")
                result.append("begin_option")
                for key, value in options.items():
                    if isinstance(value, list):
                        value = list(map(str, value))
                        result.append(" begin_{}".format(key))
                        result.append(" "+" ".join(value))
                        result.append(" end_{}".format(key))
                    else:
                        value = str(value)
                        result.append(" "+" ".join([key+"=", value]))
                result.append("end_option")
                result.append("")
            print("debug, option", result)
            return result

        card = []
        card += from_keylist(dic, ["go", "potentialfile"])
        if dic["brvtyp"] == "aux":
            card += from_keylist(dic, ["brvtyp"])
            card += make_lattice_vectors(dic)
        else:
            card += from_keylist(dic, ["brvtyp", "a",
                                       "c/a", "b/a", "alpha", "beta", "gamma"])
        card += from_keylist(dic, ["edelt", "ewidth",
                                   "reltyp", "sdftyp", "magtyp", "record"])
        card += from_keylist(dic, ["outtyp", "bzqlty", "maxitr", "pmix"])
        card += from_keylist(dic, ["ntyp"])
        card += make_type_card(dic)
        card += from_keylist(dic, ["natm"])
        card += make_atom_card(dic)
        card += make_option_card(dic)
        if dic["go"] == "fsm":
            card += [str(dic["fspin"])]
        elif dic["go"][:3] == "spc":
            if "kpath_raw":
                card += dic["kpath_raw"]
        with open(self.path_dir+"/"+inputcard, mode="w") as f:
            f.write("\n".join(card))

    def run(self, akaikkr_exe, infile, outfile):
        """run akaikkr

        stopped by errtrp is also treaed as KKRFailedExecutionError

        Args:
            akaikkr_exe (str): specx file apth
            infile (str): input filename
            outfile (str): output filename

        Raises:
            KKRFailedExecutionError: failed to execute specx or error occured in executing specx.
        Retruns:
            int: return code of akaikkr_exe
        """
        # execute akaikkr
        cmd = "cd {}; {} < {} > {}".format(
            self.path_dir, akaikkr_exe, infile, outfile)
        ret = subprocess.call(cmd, shell=True)
        if ret != 0:
            raise KKRFailedExecutionError("return_code={}".format(ret))

        stoppped_by_errtrp = self.check_stopped_by_errtrp(outfile)
        if stoppped_by_errtrp:
            raise KKRFailedExecutionError("stopped by errtrp")

        return ret

    def copy_potential(self, potentialfile1, potentialfile2,
                       path_dir1=None, path_dir2=None):
        """copy potential file

            If path_dir1 is not None, file1 =  path_dir1/potentialfile1,
            else file1 = self.directory/potentialfile1.
            If path_dir2 is not None, file2 =  path_dir2/potentialfile1,
            else file2 = self.directory/potentialfile2.

            raise errors in shutil.copyfile

        Args:
            potentialfile1 (str): input potential filename.
            potentialfile2 (str): output potential filename.
            path_dir1 (str, optional): path for potential file1. Defaults to None.
            path_dir2 (str, optional): path for potential file2. Defaults to None.
        """
        # copy potential file: file1 => file2
        if True:
            if path_dir1 is not None:
                filepath1 = os.path.join(path_dir1, potentialfile1)
            else:
                filepath1 = os.path.join(self.path_dir, potentialfile1)
            if path_dir2 is not None:
                filepath2 = os.path.join(path_dir2, potentialfile2)
            else:
                filepath2 = os.path.join(self.path_dir, potentialfile2)
            shutil.copyfile(filepath1, filepath2)
        else:
            cmd = "cd {}; cp {} {}".format(
                self.path_dir, potentialfile1, potentialfile2)
            subprocess.call(cmd, shell=True)

    def delete_potential(self, potentialfile):
        """delete potential file

            raise errors in os.remove.

        Args:
            potentialfile (str): potential file
        """
        # delete potential file
        if True:
            filepath = os.path.join(self.path_dir, potentialfile)
            os.remove(filepath)
        else:
            cmd = "cd {}; rm -rf {}".format(self.path_dir, potentialfile)
            subprocess.call(cmd, shell=True)

    def check_stopped_by_errtrp(self, outfile):
        """check whether program is stopped by errtrp

        Args:
            outfile (str): output filename

        Returns:
            bool: stoppedby errtrp
        """
        data = self._read(outfile)
        line = data[-1]
        if line.startswith(" ***err in"):
            return True
        else:
            return False

    def check_convergence_go(self, outfile):
        """check convergence for go calculation.
        converged if sbtime report is found.
        not converged if *** no convergence is found.

        Args:
            outfile (str): output filename

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            bool: converged or not
        """
        # check convergence for go calcualtion
        data = self._read(outfile)

        flag = True
        converged = False
        for line in data:
            if "*** no convergence" in line:
                flag = False
                break
            if "sbtime report" in line:
                flag = False
                converged = True
        if flag:
            raise KKRValueAquisitionError("failed to get convergence")
        return converged

    def get_threads_openmp(self, outfile):
        """get the number of thread for OpenMP

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            int: the number of threads
        """
        threads = None
        data = self._read(outfile)

        for line in data:
            if "threads" in line:
                threads = int(line.split("(")[1].split()[0])
        if threads is None:
            print("no threads in file, but continue", outfile)
        return threads

    def get_rms_error(self, outfile):
        """get rms error of iterations

        Args:
            outfile (str): output filename to read

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float]: rms errors
        """
        data = self._read(outfile)

        rms_error = []
        for line in data:
            if "te=" in line:
                rms_error.append(float(line.split("err=")[1].split()[0]))
        if len(rms_error) == 0:
            raise KKRValueAquisitionError("failed to get rms error history")
        return rms_error

    def get_err_history(self, outfile):
        """make err history from output

        Args:
            outfile (str): output of kkr

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float]: RMS errors history
        """
        return self.get_rms_error(outfile)

    def get_te_history(self, outfile):
        """get te of iterations

        Args:
            outfile (str): output filename to read

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float]: te history
        """
        # get information of rms-error
        data = self._read(outfile)

        te = []
        for line in data:
            if "te=" in line:
                te.append(float(line.split("te=")[1].split()[0]))
        if len(te) == 0:
            raise KKRValueAquisitionError("failed to get te history")
        return te

    def get_moment_history(self, outfile):
        """get moment of iterations

        Args:
            outfile (str): output filename to read

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float]: moment history
        """
        # get information of rms-error
        data = self._read(outfile)

        te = []
        for line in data:
            if "te=" in line:
                te.append(float(line.split("moment=")[1].split()[0]))
        if len(te) == 0:
            raise KKRValueAquisitionError("failed to get moment history")
        return te

    def get_lattice_constant(self, outfile, unit="au"):
        """get lattice constant (a)
        The geometries are scaled by a in Akaikkr.

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: lattice constant
        """
        # get lattice constant
        data = self._read(outfile)

        lattice_constant = None
        for line in data:
            if "bravais=" in line:
                lattice_constant = float(line.split("a=")[1].split()[0])
        if lattice_constant is None:
            raise KKRValueAquisitionError("failed to get lattice_constant")
        if unit.lower() == "ang":
            lattice_constant *= _au2ang
        return lattice_constant

    def get_struc_param(self, outfile, unit="au"):
        """get lattice constant (a)
        The geometries are scaled by a in Akaikkr.

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: lattice constant
        """
        # get lattice constant
        data = self._read(outfile)

        lattice_constant = None
        coa = None
        boa = None
        alpha = None
        beta = None
        gamm_ = None
        for line in data:
            if "bravais=" in line:
                brvtyp = line.split("bravais=")[1].split()[0]
                lattice_constant = float(line.split("a=")[1].split()[0])
                coa = float(line.split("c/a=")[1].split()[0])
                boa = float(line.split("b/a=")[1].split()[0])
            if "alpha=" in line:
                alpha = float(line.split("alpha=")[1].split()[0])
                beta = float(line.split("beta=")[1].split()[0])
                gamma_ = float(line.split("gamma=")[1].split()[0])
        if lattice_constant is None or \
           brvtyp is None or coa is None or boa is None or \
           alpha is None or beta is None or gamma_ is None:
            raise KKRValueAquisitionError("failed to get struc_param")
        if unit.lower() == "ang":
            lattice_constant *= _au2ang
            coa *= _au2ang
            coa *= _au2ang

        return {"brvtyp": brvtyp, "a": lattice_constant, "c/a": coa,
                "b/a": boa, "alpha": alpha, "beta": beta, "gamma": gamma_}

    def get_ntype(self, outfile):
        """get ntype in the output file

        Args:
            outfile (str): an akaikkr output filename

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            int: ntype
        """
        data = self._read(outfile)
        ntyp = None
        for line in data:
            if line.startswith("   ntyp="):
                ntyp = int(line.split("ntyp=")[1].split()[0].strip())
                break
        if ntyp is None:
            raise KKRValueAquisitionError("failed to get ntyp")
        return ntyp

    if False:
        # module related to dsplc isn't made of output.
        # They are placed in DsplcMaker class.
        def get_dsplc(self, structure, alen, cart_len=0.01, direction=None,
                      frac=True):
            """get a list of displc.

            displc is scaled by alen.

            Args:
                structure (Structure): pymatgen.core.Structure
                alen (float): a length
                cart_len (float, optional): lenght of dx in cartesian. Defaults to 0.01.
                direction (nd.ndarray, optional): initial displacement vector. Defaults to None.

            Returns:
                [[float,float,float]]: rotated/operated cartesian coordinates
            """
            dsplcmaker = DsplcMaker(structure)
            displc = dsplcmaker.get_displc_as_cartesian(
                cart_len=cart_len, direction=direction, frac=frac)

            # scale displc by alen
            displc = np.array(displc)/alen
            displc = displc.tolist()

            return displc

        def get_avr_debye(self, m, debye_temp, temp):
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

    def conversion_block_list(self, block_list):
        """conversion helper routine from block_list to anclrXX.concYY

        Args:
            block_list (list): a block list of types

        Returns:
            [Str]: type list
        """
        def add_component_name(complist):
            """add "comp_shortname" to comlist dict

            Args:
                complist (list): a list of dict

            Returns:
                list: a list of dict where "comp_shortname" elements are added.
            """
            for comp in complist:
                name = []
                if len(comp["component"]) > 1:  # alloy
                    for x in comp["component"]:
                        Z = int(x["anclr"])
                        try:
                            symbol = Element("H").from_Z(Z)
                            s = "{}_{}_{:.1f}%".format(
                                comp["type"], symbol, 100*x["conc"])
                        except ValueError:
                            s = "{}_Z{}_{:.1f}%".format(
                                comp["type"], Z, 100*x["conc"])
                        name.append(s)
                else:  # stoichiometry
                    name = [comp["type"]]
                comp["comp_shortname"] = name
            return complist

        type_list = []
        for block in block_list:
            block_dic = {"component": None}
            for line in block:
                if "type=" in line:
                    typename = line.split("type=")[1].rstrip().split()[0]
                    block_dic["type"] = typename
                    lmx = line.split("lmxtyp=")[1].rstrip().split()[0]
                    block_dic["lmx"] = lmx
                elif "anclr=" in line:
                    # component_id = line.split("cmpnt")[1].rstrip().split()[0]
                    anclr = line.split("anclr=")[1].rstrip().split()[0]
                    # read x= or conc=
                    if "x=" in line:
                        conc = line.split("x=")[1].rstrip().split()[0]
                    elif "conc=" in line:
                        conc = line.split("conc=")[1].rstrip().split()[0]
                    component_dic = {"anclr": float(
                        anclr), "conc": float(conc)}
                    if block_dic["component"] is None:
                        block_dic["component"] = []
                    block_dic["component"].append(component_dic)
            type_list.append(block_dic)
        if len(type_list) == 0:
            raise KKRValueAquisitionError(
                "failed to get type of site, conversion_block_list")
        type_list = add_component_name(type_list)
        return type_list

    def get_type_of_site(self, outfile, ntype=None):
        """get the type of site part. automatically read ntype if ntype is None.

        Args:
            outfile (str):  an akaikkr output filename
            ntype (int, optional): number of types. Defaults to None.

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [str]: a list of types
        """
        def get_block_list(data):
            """ get clock_list from data
            """
            istart = None
            for i, line in enumerate(data):
                if line.startswith("   type of site"):
                    istart = i
                    break
            for i in range(istart+1, len(data)):
                line = data[i].rstrip()
                if len(line) == 0:
                    iend = i
                    break
            # data[istart+1:iend]
            block_list = []
            block = []
            for line in data[istart+1:iend]:
                if line.startswith("   type="):
                    if len(block) > 0:
                        block_list.append(block)
                        block = []
                block.append(line.rstrip())
            if len(block) > 0:
                block_list.append(block)
            return block_list

        if ntype is None:
            ntype = self.get_ntype(outfile)
        filepath = os.path.join(self.path_dir, outfile)
        with open(filepath) as f:
            data = f.readlines()
        block_list = get_block_list(data)
        if len(block_list) == 0:
            raise KKRValueAquisitionError("failed to get type of site")
        return self.conversion_block_list(block_list)

    def cut_pdos_all(self, dosfile, output_format="sequential"):
        """cut PDOS

          '-2.0550       0.0219    0.0598    0.089912188.5942' occurs.

        Thus, first floating point format is analyzed.

        output_format = "sequential|spin_separation"
        output_format=sequential output pdos as the order of outputcard
        output_format=spin_seperation output pdos, where the first index is for spin up/down

        Args:
            dosfile (str): an akaikkr output filename
            output_format (bool, optional): output format. Defaults to "sequntial"

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float], [[[float]]]: energy mesh, up/down spin PDOS
        """
        data = self._read(dosfile)
        keyword = "DOS of component"

        value_positions_list = _analyze_pdos_format(keyword, data)
        energy = None

        if True:
            pdos_up_list = []
            pdos_dn_list = []
            for icmp, value_positions in enumerate(value_positions_list):
                data_iter = iter(data)
                search_keyword = " {}{:2d}".format(keyword, icmp+1)

                e, pdos = _cut_a_pdos_block(
                    data_iter, search_keyword, value_positions)
                pdos_up_list.append(pdos)
                if e is not None:
                    energy = e
                e, pdos = _cut_a_pdos_block(
                    data_iter, search_keyword, value_positions)
                if pdos is not None:
                    pdos_dn_list.append(pdos)

            if len(pdos_up_list) == 0:
                raise KKRValueAquisitionError("failed to get pdos")

            if output_format == "spin_separation":
                if len(pdos_dn_list) > 0:
                    pdos_block = [pdos_up_list, pdos_dn_list]
                else:
                    pdos_block = [pdos_up_list]
            elif output_format == "sequential":
                pdos_block = []
                for pdos in pdos_up_list:
                    pdos_block.append(pdos)
                for pdos in pdos_dn_list:
                    pdos_block.append(pdos)

            if len(pdos_block) == 0:
                raise KKRValueAquisitionError("failed to get pdos")

            return energy, pdos_block
        else:
            icmp = 0
            pdos_block = []
            e = []
            pdos = []
            check = False
            for line in data:
                if keyword in line:
                    check = True
                    continue
                if check and len(line.strip()) == 0:
                    check = False
                    pdos_block.append(pdos)
                    energy = e
                    e = []
                    pdos = []
                    icmp += 1
                    continue
                if check:
                    value_positions = value_positions_list[icmp]
                    values = [line[istart:iend]
                              for istart, iend in value_positions]
                    e.append(float(values[0]))
                    v2 = []
                    for s in values[1:]:
                        if s.startswith("**********"):
                            v2.append("1000.0")
                        else:
                            v2.append(s)
                    try:
                        v = list(map(float, v2))
                    except ValueError:
                        print("\nerror in cutpdos_all, path_dir=",
                              self.path_dir, "filename=", dosfile)
                        print("return None, None\n")
                        raise KKRValueAquisitionError("failed to get pdos")
                        return None, None
                    pdos.append(v)
            if energy is None:
                print("get_pdos_all() failed. path=", self.path_dir, dosfile)
                print("return None, None\n")
                raise KKRValueAquisitionError("failed to get pdos")
                return None, None
        return energy, pdos_block

    def make_pdos_all_df(self, dosfile):
        """get pdos all as DataFrame
        output contains dict["up"], and dict["dn"]

        Args:
            dosfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            dict:  a dictionary of DataFrame
        """
        def pdos2df(energy, up):
            df_energy = pd.DataFrame(energy, columns=["energy"])
            up = np.array(up)
            df_up = pd.DataFrame(up, columns=l_label[:up.shape[1]])
            df_pdosup = pd.concat([df_energy, df_up], axis=1)
            return df_pdosup

        updn_label = ["up", "dn"]  # shoule be outside of the member.
        l_label = ["s", "p", "d", "f", "g"]
        outputformat = "spin_separation"
        energy, pdos_block = self.cut_pdos_all(
            dosfile, output_format=outputformat)
        pdos_up_list = []
        pdos_dn_list = []
        if len(pdos_block) > 1:
            for up, dn in zip(pdos_block[0], pdos_block[1]):
                df_pdosup = pdos2df(energy, up)
                pdos_up_list.append(df_pdosup)
                df_pdosdn = pdos2df(energy, dn)
                pdos_dn_list.append(df_pdosdn)
            return {updn_label[0]: pdos_up_list, updn_label[1]: pdos_dn_list}
        else:
            for up in pdos_block[0]:
                df_pdosup = pdos2df(energy, up)
                pdos_up_list.append(df_pdosup)
            return {updn_label[0]: pdos_up_list, updn_label[1]: []}

    def get_magtyp(self, outfile):
        """get magtype

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: magtype
        """
        data = self._read(outfile)
        magtype = None
        for line in data:
            if "magtyp=" in line:
                magtype = line.split("magtyp=")[1].split()[0]
        if magtype is None:
            raise KKRValueAquisitionError("failed to get magtyp")
        return magtype

    def get_unitcell_volume(self, outfile):
        """get the unit cell volume

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: the unit cell volume
        """
        # get volume
        data = self._read(outfile)

        unitcell_volume = None
        for line in data:
            if "unit cell volume=" in line:
                unitcell_volume = float(line.rstrip(
                    "\n").split("volume=")[1][:-6])
        if unitcell_volume is None:
            raise KKRValueAquisitionError("failed to get unit cell volume")
        return unitcell_volume

    def get_ewidth(self, outfile):
        """get ewidth, which is the minimum values of energy integration.

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: ewidth
        """
        # get ewidth
        data = self._read(outfile)
        ewidth = None
        for line in data:
            if "ewidth=" in line:
                ewidth = float(line.split("ewidth=")[1].split()[0])
        if ewidth is None:
            raise KKRValueAquisitionError("failed to get ewidth")
        return ewidth

    def get_go(self, outfile):
        """get go keyword

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            str: go keyword
        """
        data = self._read(outfile)
        go = None
        for line in data:
            if "go=" in line:
                go = line.split("go=")[1].split()[0]
        if go is None:
            raise KKRValueAquisitionError("failed to get go")
        return go

    def get_potentialfile(self, outfile):
        """get potentialfile

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            str: potentialfile
        """
        data = self._read(outfile)
        f = None
        for line in data:
            if "file=" in line:
                f = line.split("file=")[1].split()[0]
        if f is None:
            raise KKRValueAquisitionError("failed to get potential file")
        return f

    def get_edelt(self, outfile):
        """get edelt, imaginary part of complex integration path

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: edelt
        """
        # get edelt
        data = self._read(outfile)
        edelt = None
        for line in data:
            if "edelt=" in line:
                edelt = float(line.split("edelt=")[1].split()[0])
        if edelt is None:
            raise KKRValueAquisitionError("failed to get edelt")
        return edelt

    def get_fermi_level(self, outfile):
        """get fermi level.
            the average in the case of mag.

        Args:
            outfile (str): output filename

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: fermi level
        """
        # get Fermi level
        data = self._read(outfile)
        ef = None
        for line in data:
            if "ef=" in line:
                ef_up = float(line.split()[1])
                ef_dn = float(line.split()[2])
                ef = (ef_up+ef_dn)/2
                break
        if ef is None:
            raise KKRValueAquisitionError("failed to get ef")
        return ef

    def get_total_energy(self, outfile):
        """get total energy (te)

        Args:
            outfile (st): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: total energy (te)
        """
        # get total energy
        data = self._read(outfile)
        te = None
        for line in data:
            if "total energy=" in line:
                te = float(line.split("total energy=")[-1].rstrip())
                break
        if te is None:
            raise KKRValueAquisitionError("failed to get te")
        return te

    def get_fixed_spin_moment(self, outfile):
        """get fixed spin moment

        Args:
            outfile (str): output filename

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: fixed spin moment
        """
        # get total energy
        data = self._read(outfile)
        fspin = None
        for line in data:
            if "fixed spin moment=" in line:
                fspin = float(line.split("fixed spin moment=")[1].rstrip())
                break
        if fspin is None:
            raise KKRValueAquisitionError("failed to get fspin")
        return fspin

    def get_total_moment(self, outfile):
        """get total moment

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: the final total moment
        """
        # get total moment
        data = self._read(outfile)
        last_moment = None
        h_moment = []
        for line in data:
            if "te=" in line:
                h_moment.append(float(line.split("moment=")[1].split()[0]))
        if len(h_moment) > 0:
            last_moment = h_moment[-1]
        else:
            raise KKRValueAquisitionError("failed to get total moment")
        return last_moment

    def get_local_moment(self, outfile, mode="spin"):
        """get local moment.
        There are spin and orbital contributions.
        output is spin moments and orbital moments if mode=="all".
        output is spin moments if mode=="spin".
        output is orbital moments if mode=="orbital".

        Args:
            outfile (str): filename to analyze
            mode (str, optional): output mode. Defaults to "spin".

        Raises:
            ValueError: mode is unknown
            KKRValueAquisitionError: failed to get keyword

        Returns:
            list: a list of spin moment value
        """
        # get information of local moment(spin or orbital)
        data = self._read(outfile)

        spin_moment = []
        orbital_moment = []
        for line in data:
            if "orbital moment=" in line:
                spin_moment.append(
                    float(line.split("spin moment=")[1].split()[0]))
                orbital_moment.append(
                    float(line.split("orbital moment=")[1].split()[0]))
        if len(spin_moment) == 0:
            raise KKRValueAquisitionError("failed to get spin moment")
            spin_moment = None
        if len(orbital_moment) == 0:
            raise KKRValueAquisitionError("failed to get orbital moment")
            orbital_moment = None
        if mode == "spin":
            return spin_moment
        elif mode == "orbital":
            return orbital_moment
        elif mode == "all":
            return spin_moment, orbital_moment
        else:
            raise ValueError("failed to get local moment. "
                             "because of unknown mode={}".format(mode))

    def get_type_charge(self, outfile):
        """get type charges

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            list: charge values of types
        """
        # get total charge of type
        data = self._read(outfile)

        charge = []
        for line in data:
            if "total charge=" in line:
                charge.append(float(line.split("total charge=")[1].split()[0]))
        if len(charge) == 0:
            raise KKRValueAquisitionError("failed to get charge")
            charge = None
        return charge

    def get_prim_vec(self, outfile, unitof="relative"):
        """get primitive vector.
        output is as an unit of a-length if unitof=="relative".

        Args:
            outfile (str): filename to analyze
            unitof (str, optional): the unit of primitive vectors. Defaults to "relative".

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            list: primitive vectors
        """
        data = self._read(outfile)

        data_iter = iter(data)
        for line in data_iter:
            line = line.strip()
            s = line.split()
            if len(s) > 3:
                if s[0] == "primitive" and s[1] == "translation" and s[2] == "vectors":
                    break
        prim_vec = []
        for i in range(3):
            line = next(data_iter)
            line = line.replace("a=(", "").replace(
                "b=(", "").replace("c=(", "").replace(")", "")
            s = line.strip().split()
            v = list(map(float, s))
            prim_vec.append(v)
        if len(prim_vec) == 0:
            raise KKRValueAquisitionError("failed to get primitive vector")
        return prim_vec

    def get_atom_coord(self, outfile,):
        """get atomic coordinates and names

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            list, list: atom coodinates, atom names
        """
        data = self._read(outfile)

        data_iter = iter(data)
        for line in data_iter:
            line = line.strip()
            s = line.split()
            if len(s) > 3:
                if s[0] == "atoms" and s[1] == "in" and s[2] == "the":
                    break
        atom_coords = []
        atom_names = []
        while True:
            line = next(data_iter)
            s = line.strip().split()
            if len(s) == 0:
                break
            if s[0] == "position=":
                v = list(map(float, s[1:4]))
                atom_coords.append(v)
                name = "".join(s[4:]).replace("type=", "")
                atom_names.append(name)
            else:
                break
        if len(atom_coords) == 0:
            raise KKRValueAquisitionError("failed to get atom_coords")
            atom_coords = None
        if len(atom_names) == 0:
            raise KKRValueAquisitionError("failed to get atom_names")
            atom_names = None
        return atom_coords, atom_names

    def make_pymatgenstructure(self, outfile, unit="ang", change_atom_name=False):
        """make pymatgen Structure from output

        Args:
            outfile (str): outputcard filename
            unit (str, optional): unit of length. Defaults to "ang".
            change_atom_name (bool, optional): change atomic name to Element names. Defaults to False.

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            pymatgen.Structure: structure

        """
        def _find_component(name, typeofsites):
            for _site in typeofsites:
                if name == _site["type"]:
                    return _site["component"]
            print("failed to find type", name)
            raise KKRValueAquisitionError(
                "Error in make_pymatgenstructure, failed to find type")

        alen = self.get_lattice_constant(outfile, unit=unit)
        prim_vec = self.get_prim_vec(outfile)
        prim_vec = np.array(prim_vec)*alen
        atom_coords, atom_names = self.get_atom_coord(outfile)
        atom_coords = np.array(atom_coords)*alen
        typeofsites = self.get_type_of_site(outfile)

        if change_atom_name:
            sites = _noalloy_atom_names(atom_names)
        else:
            sites = []
            for name in atom_names:
                component = _find_component(name, typeofsites)
                pymatgensite = {}
                for comp in component:
                    anclr = comp["anclr"]
                    Z = int(anclr)
                    if Z > 0:
                        elm = str(Element("H").from_Z(Z))
                        occu = comp["conc"]
                        pymatgensite.update({elm: occu})

                sites.append(pymatgensite)

        struc = Structure(lattice=prim_vec, species=sites, coords=atom_coords,
                          coords_are_cartesian=True)
        return struc

    def make_poscar(self, outfile, change_atom_name=False, poscar_filename=None):
        """make poscar.
        make a file if poscar_filename is not None.

        This function will become obsolete.

        Args:
            outfile (str): filename to analyze
            change_atom_name (str, optional):: change atom name to Elements. Defaults to False
            poscar_filename (str, optional): poscar filename. Defaults to None.

        Returns:
            [str]: poscar as a list of string
        """
        _msg = """This function will become obsolete.
Please consider using
struc = job.make_pymatgenstructure(outfile)
StructureWriter.poscar(struc)."""
        print(_msg)

        alen = self.get_lattice_constant(outfile)
        prim_vec = self.get_prim_vec(outfile)
        atom_coords, atom_names = self.get_atom_coord(outfile)

        uniq_atom_names = list(set(atom_names))
        name_coords_list = []
        for uniq_name in uniq_atom_names:
            name_coords = []
            for coords, name in zip(atom_coords, atom_names):
                if name == uniq_name:
                    name_coords.append(coords)
            name_coords_list.append(name_coords)

        # pymatgen can't read unknown name
        # Thus element name must be corrected.
        if change_atom_name:
            uniq_atom_names = _fix_uniq_atom_names(uniq_atom_names)

        # make poscar outputs
        s = []
        for uniq_name, name_coords in zip(uniq_atom_names, name_coords_list):
            n = len(name_coords)
            s.append("{}{}".format(uniq_name, n))
        compound_name = "".join(s)
        lines = [compound_name]
        lines.append("{}".format(alen*0.529177))
        for v in prim_vec:
            s = list(map(str, v))
            lines.append(" ".join(s))
        lines.append(" ".join(uniq_atom_names))
        iv = [len(v) for v in name_coords_list]
        s = list(map(str, iv))
        lines.append(" ".join(s))
        lines.append("cartesian")
        for uniq_name, name_coords in zip(uniq_atom_names, name_coords_list):
            for coords in name_coords:
                s = list(map(str, coords))
                lines.append(" ".join(s))
        if poscar_filename is not None:
            with open(poscar_filename, "w") as f:
                f.write("\n".join(lines))
            # print("  saved to", poscar_filename)
        return lines

    def get_curie_temperature(self, outfile):
        """get curie temperature

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: curie temperature
        """
        # get curie temperature
        data = self._read(outfile)
        tc = None
        for line in data:
            if "Tc (in mean field approximation) =" in line:
                tc = float(line.split(
                    "Tc (in mean field approximation) =")[-1].rstrip()[:-1].strip())
                break
        if tc is None:
            raise KKRValueAquisitionError("failed to get curie temperature")
        return tc

    def get_resistivity(self, outfile):
        """get resistivity

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            float: resistivity
        """
        # get resistivity
        data = self._read(outfile)
        resistivity = None
        for line in data:
            if "resistivity" in line:
                resistivity = float(line.split()[1])
                break
        if resistivity is None:
            raise KKRValueAquisitionError("failed to get resistivity")
        return resistivity

    def get_conductivity_spin(self, outfile):
        """get conductivity per spin

        Args:
            outfile (str): filename to analyze

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            list: a list of conductivity per spin
        """
        # get resistivity
        data = self._read(outfile)

        conductivity_spin = []
        for line in data:
            if "cnd1" in line:
                conductivity_spin.append(float(line.split()[1]))
            if "cnd2" in line:
                conductivity_spin.append(float(line.split()[1]))
                break
        if len(conductivity_spin) == 0:
            raise KKRValueAquisitionError("failed to get cnd")
            conductivity_spin = None
        return conductivity_spin

    def check_core_level(self, outfile, core_state=["3d", "4d", "4f"]):
        """check core states defined by core_state.

        Args:
            outfile (str): filename to analyze
            core_state ([str]): core states

        Returns:
            list: a list where the core exist or not
            list: a list of the core levels
        """
        # check core states
        core_exist = []
        core_level = []
        data = self._read(outfile)

        for core in core_state:
            dummy1 = False
            dummy2 = []
            for line in data:
                if "Ry("+core+")" in line:
                    dummy1 = True
                    dummy2.append(
                        float(line.split("Ry("+core+")")[0].split()[-1]))
            core_exist.append(dummy1)
            core_level.append(dummy2)

        # accept that they are none
        return core_exist, core_level

    def cut_dos_float(self, dosfile, keyword="total DOS", save=False):
        """get DOS as float values.
        The filename of the DOS csv file is keyword.replace(" ","_").csv.

        Args:
            dosfile (str): filename to analyze
            keyword (str, optional): keyword to find DOS. Defaults to "total DOS".
            save (bool, optional): save DOS file or not as csv. Defaults to False.

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float],[float],[float]: energies, DOS up, DOS down
        """
        # cut dos (beta version)
        data = self._read(dosfile)

        e = []
        dos_up = []
        dos_dn = []
        nspin = 0
        check = False
        for line in data:
            line = line.strip()
            if keyword in line:
                check = True
                nspin += 1
            if not line:
                check = False
            if check and keyword not in line:
                if nspin == 1:
                    e.append(float(line.split()[0]))
                    dos_up.append(float(line.split()[1]))
                else:
                    dos_dn.append(float(line.split()[1]))

        if len(e) == 0 or len(dos_up) == 0:
            raise KKRValueAquisitionError("failed to get dos")

        if save:
            filename = os.path.join(
                self.path_dir, keyword.replace(" ", "_")+".csv")
            if nspin == 2:
                df = pd.DataFrame(
                    {"energy": e, "dos_up": dos_up, "dos_dn": dos_dn})
            else:
                df = pd.DataFrame(
                    {"energy": e, "dos_up": dos_up, })
            df.to_csv(filename, index=False)
            print("save to", filename)

        return e, dos_up, dos_dn

    def cut_dos(self, dosfile, keyword="total DOS", save=False):
        """get DOS as float values (same as cut_dos_float()).
        The filename of the DOS csv file is keyword.replace(" ","_").csv.

        Args:
            dosfile (str): filename to analyze
            keyword (str, optional): keyword to find DOS. Defaults to "total DOS".
            save (bool, optional): save DOS file or not as csv. Defaults to False.

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [float],[float],[float]: energies, DOS up, DOS down
        """
        self.cut_dos_float(dosfile, keyword, save)

    def cut_dos_raw(self, dosfile, keyword, save=False):
        """get DOS as lines.
        Please use this routine if cut_dos_float() doesn't work.

        Args:
            dosfile (str): filename to analyze
            keyword (str): keyword to find DOS
            save (bool, optional): save the DOS csv file or not. Defaults to False.

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            [str], [str]: a list of DOS up lines, a list of DOS down lines
        """
        # cut dos
        data = self._read(dosfile)

        dos_up = []
        dos_dn = []
        nspin = 0
        check = False
        for line in data:
            line = line.strip()
            if keyword in line:
                check = True
                nspin += 1
            if not line:
                check = False
            if check and keyword not in line:
                if nspin == 1:
                    dos_up.append(line)
                else:
                    dos_dn.append(line)

        if (dos_up) == 0:
            raise KKRValueAquisitionError("failed to get dos in the raw mode")

        if save:
            with open(self.path_dir+"/"+keyword.replace(" ", "_")+"_up.dat", "w") as f:
                f.write("\n".join(dos_up)+"\n")
            if len(dos_dn) > 0:
                with open(self.path_dir+"/"+keyword.replace(" ", "_")+"_dn.dat", "w") as f:
                    f.write("\n".join(dos_dn)+"\n")
        else:
            return dos_up, dos_dn

    def _jij_lines_dataframe(self, lines):
        """convert jij lines to dataframe

        Args:
            lines ([str]): jij section of kkr outputcard

        Returns:
            pd.DataFrame: jij
        """
        label = lines[0][0]
        columns = lines[1]

        i = columns.index("site")
        del columns[i]
        columns.insert(i, "site1")
        columns.insert(i+1, "site2")
        i = columns.index("comp")
        del columns[i]
        columns.insert(i, "comp1")
        columns.insert(i+1, "comp2")
        i = columns.index("cell")
        del columns[i]
        columns.insert(i, "a")
        columns.insert(i+1, "b")
        columns.insert(i+2, "c")

        df = pd.DataFrame(lines[2:], columns=columns)
        labels = label.split("-")
        label1, label2 = labels[0], labels[1]
        df["type1"] = label1
        df["type2"] = label2
        return df

    def cut_jij_dataframe(self, jijfile):
        """cut jij section and make dataframe

        Args:
            jijfile (str): outputfile name

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            pd.DataFrame: jij
        """
        data = self._read(jijfile)

        data_iter = iter(data)

        # find J_ij
        while True:
            line = next(data_iter)
            s = line.strip()
            if s == "J_ij":
                break
        # split by empty lnnes
        all_lines = []
        lines = []
        while True:
            line = next(data_iter).strip()
            if line.startswith("Tc"):
                break
            elif len(line) == 0:
                all_lines.append(lines)
                lines = []
            else:
                lines.append(line.split())
        if len(lines) > 0:
            raise KKRValueAquisitionError("failed to get jij")

        df_list = []
        for lines in all_lines:
            df = self._jij_lines_dataframe(lines)
            df_list.append(df)
        df = pd.concat(df_list, axis=0)
        del df["index"]  # delete the original index of the file
        return df

    def cut_jij_raw(self, jijfile, keyword, comp1=1, comp2=1, save=False):
        """get Jij.
        The filename of the Jij is keyword.replace(" ","_").dat.

        Args:
            jijfile (str): filename to analyze
            keyword (str): keyword to find Jij part in the filename
            comp1 (int, optional): comp1 index. Defaults to 1.
            comp2 (int, optional): comp2 index. Defaults to 1.
            save (bool, optional): save extracted data or not. Defaults to False.

        Raises:
            KKRValueAquisitionError: failed to get keyword

        Returns:
            list: a list of string lines only for (comp1,comp2)
        """
        # cut jij
        data = self._read(jijfile)

        jij = []
        check = False
        for line in data:
            line = line.strip()
            if keyword in line:
                check = True
            if not line:
                check = False
            if check and keyword not in line and "index" not in line:
                if int(line.split()[3]) == comp1 and int(line.split()[4]) == comp2:
                    jij.append(line)
        if len(jij) == 0:
            raise KKRValueAquisitionError("failed to get jij")
        if save:
            outfile = self.path_dir+"/"+keyword + \
                "_"+str(comp1)+"-"+str(comp2)+".dat"
            with open(outfile, "w") as f:
                f.write("\n".join(jij)+"\n")
        else:
            return jij

    def show_summary_go(self, outfile):
        """show human readable summary

        Args:
            outfile (str): output filename of akaikkr to analyze
        """
        go = self.get_go(outfile)

        # get summary of go calculation
        print("  directory: ", self.path_dir, "  output file: ", outfile)
        print("  converged: ", self.check_convergence_go(outfile))
        print("  threads: ", self.get_threads_openmp(outfile))
        if go == "fsm":
            print("  fixed spin moment: ",
                  self.get_fixed_spin_moment(outfile), " muB")
        print("  total energy: ", self.get_total_energy(outfile), " Ry")
        print("  total moment: ", self.get_total_moment(outfile), " muB")
        print("  local moment: ", self.get_local_moment(outfile), " muB")
        if go == "tc" or go[0] == "j":
            print("  Tc: ", self.get_curie_temperature(outfile), "K")
        if go == "dos":
            ene, dos_up, dos_dn = self.cut_dos_float(outfile)
            ene = np.array(ene)
            dos_up = np.array(dos_up)
            dos_dn = np.array(dos_dn)
            idx = np.abs(ene).argsort()
            iv = idx[:2]
            if len(dos_dn) > 0:
                for i in iv:
                    print("  up/dn total DOS (energy):  {} {} ({} Ry)".format(
                        dos_up[i], dos_dn[i], ene[i]))
            else:
                for i in iv:
                    print("  total DOS (energy):  {} ({} Ry)".format(
                        dos_up[i], ene[i]))
        if go == "cnd":
            print("  Resistivity: ", self.get_resistivity(
                outfile), " micro ohm cm")
        print()

    def get_result_testrun(self, outfile, Awk_k=None):
        """get result for testrun

        Args:
            outfile (str): output filename of akaikkr to analyze
            Awk_k (int, optional): k index to save A(w,k). Defaults to None.

        Returns:
            dict: result summary
        """

        go = self.get_go(outfile)

        dic = {"go": go, "te": self.get_total_energy(outfile),
               "rms": self.get_rms_error(outfile)[-1],
               "tm": self.get_total_moment(outfile),
               "spinlocalmoment": self.get_local_moment(outfile, mode="spin"),
               "orbitallocalmoment": self.get_local_moment(outfile, mode="orbital"),
               "threads": self.get_threads_openmp(outfile),
               "conv": self.check_convergence_go(outfile)}
        if go == "fsm":
            dic.update({"fspin": self.get_fixed_spin_moment(outfile)})
        if go[0] == "j":
            dic.update({"Tc": self.get_curie_temperature(outfile)})
            _jijdf = self.cut_jij_dataframe(outfile)
            jijselect = _cut_only_nn(_jijdf)
            dic.update({"jij": jijselect})
        if go == "tc":
            dic.update({"Tc": self.get_curie_temperature(outfile)})
        if go == "dos":
            ene, dos_up, dos_dn = self.cut_dos_float(outfile)
            dos_dic = {"energy": ene, "up": dos_up, "dn": dos_dn}
            dic.update({"totalDOS": dos_dic})

            pdos = self.make_pdos_all_df(outfile)
            iatom = 0
            l_label = "d"
            updn = "up"
            df = pdos[updn][iatom]
            ene = df.loc[:, "energy"].values.tolist()
            pdos_up_d = df.loc[:, l_label].values.tolist()
            updn = "dn"
            pdos_dn_d = []
            if len(pdos[updn]) >= 1:
                df = pdos[updn][iatom]
                pdos_dn_d = df.loc[:, l_label].values.tolist()
            pdos_dic = {"energy": ene, "up": pdos_up_d, "dn": pdos_dn_d,
                        "atom": iatom, "l": l_label}
            dic.update({"pdos": pdos_dic})

        if go == "cnd":
            resistivity = self.get_resistivity(outfile)
            conduct = self.get_conductivity_spin(outfile)
            dic.update({"resis": resistivity, "cnd": conduct})
        if go[:3] == "spc":
            magtyp = self.get_magtyp(outfile)
            pot = self.get_potentialfile(outfile)
            Awk_dic = {}
            if magtyp[0] == "m":  # nspin = 2
                filename = "{}_up.spc".format(pot)
                filepath = os.path.join(self.path_dir, filename)
                Awk = AwkReader(filepath)
                if Awk_k is None:
                    Awk_dic["kpath"] = Awk.kcrt.tolist()
                    Awk_dic["kdist"] = Awk.kdist.tolist()
                    Awk_dic["energy"] = Awk.energy.tolist()
                    Awk_dic["Awk_up"] = Awk.Awk.tolist()
                    filename = "{}_dn.spc".format(pot)
                    filepath = os.path.join(self.path_dir, filename)
                    Awk = AwkReader(filepath)
                    Awk_dic["Awk_dn"] = Awk.Awk.tolist()
                else:
                    Awk_dic["energy"] = Awk.energy.tolist()
                    Awk_dic["Aw_up"] = Awk.Awk[Awk_k].tolist()
                    Awk_dic["Awk_k"] = Awk_k
                    if Awk.kpath is not None:
                        Awk_dic["Awk_kpoint"] = Awk.kpath[Awk_k]
                    else:
                        Awk_dic["Awk_kpoint"] = None
                    filename = "{}_dn.spc".format(pot)
                    filepath = os.path.join(self.path_dir, filename)
                    Awk = AwkReader(filepath)
                    Awk_dic["Aw_dn"] = Awk.Awk[Awk_k].tolist()
            else:  # nspin = 1
                filename = "{}_up.spc".format(pot)
                filepath = os.path.join(self.path_dir, filename)
                Awk = AwkReader(filepath)
                if Awk_k is None:
                    Awk_dic["kpath"] = Awk.kcrt.tolist()
                    Awk_dic["kdist"] = Awk.kdist.tolist()
                    Awk_dic["energy"] = Awk.energy.tolist()
                    Awk_dic["Awk_up"] = Awk.Awk.tolist()
                else:
                    Awk_dic["energy"] = Awk.energy.tolist()
                    Awk_dic["Aw_up"] = Awk.Awk[Awk_k].tolist()
                    Awk_dic["Awk_k"] = Awk_k
                    if Awk.kpath is not None:
                        Awk_dic["Awk_kpoint"] = Awk.kpath[Awk_k]
                    else:
                        Awk_dic["Awk_kpoint"] = None
            dic.update({"Awk": Awk_dic})
        return dic
