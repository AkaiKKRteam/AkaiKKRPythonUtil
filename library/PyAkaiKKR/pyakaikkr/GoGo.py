# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


import os
from abc import abstractmethod
import pandas as pd
import numpy as np
from pyakaikkr import plotband, KLABEL_FILENAME_, AkaikkrJob, HighSymKPath, plot_dos, plot_pdos_all, Fmg, AwkReader


def _make_displc(anclr, displc=[0, 0, 0]):
    displc_list = []
    for anclr1 in anclr:
        displc1_list = []
        for z in anclr1:
            displc1_list.append(displc)
        displc_list.append(displc1_list)
    return displc_list


def _sort_types_inside_row(df):
    """exchange type1 and type2 as type1<type2 to compare with another data.
    and add pair column.
    """
    df = df.copy().reset_index(drop=True)
    type1 = df["type1"].values
    type2 = df["type2"].values
    t1t2_list = []
    for t1, t2 in zip(type1, type2):
        _t1t2 = [t1, t2]
        _t1t2.sort()
        t1t2_list.append(_t1t2)
    df[["type1", "type2"]] = t1t2_list

    # add pair column
    comp1 = df["comp1"].values
    comp2 = df["comp2"].values
    type1 = df["type1"].values
    type2 = df["type2"].values
    typepair = []
    for t1, t2, c1, c2 in zip(type1, type2, comp1, comp2):
        typepair.append("-".join([t1, t2, c1, c2]))
    df_pair = pd.DataFrame({"pair": typepair})

    jijdf = pd.concat([df, df_pair], axis=1)
    return jijdf


def _cut_only_nn(_jijdf):
    """cut only the first columns if jij DataFrame

    Args:
        _jijdf (pd.DataFrame): Jij DataFrame

    Returns:
        pd.DataFrame: Jij only the first n.n. neighbor
    """

    jijdf = _sort_types_inside_row(_jijdf)

    typepair = np.unique(jijdf["pair"].values).tolist()
    typepair.sort()
    df_list = []
    for pair in typepair:
        _df = jijdf.query("pair=='{}'".format(pair)).sort_values(
            by="distance").reset_index(drop=True)
        _dist = _df.loc[0, "distance"]
        df = _df[_df["distance"] == _dist]
        df_list.append(df)
    dfsmall = pd.concat(df_list, axis=0)

    # in the order of distance
    dfsmall.sort_values(by="distance", inplace=True)
    dfsmall.reset_index(drop=True, inplace=True)

    jijnn = [list(dfsmall.columns)]
    v = dfsmall.values.tolist()
    jijnn.extend(v)
    return jijnn


def _make_inputcard(go):
    return "inputcard_{}".format(go)


def _make_outputcard(go):
    return "out_{}.log".format(go)


class GoGo:
    """run Akaikkr as go="go"
    """

    def __init__(self, prog,
                 directory, common_param, comment,
                 go="go", inputcard=None, outputcard=None, **args):
        """initial routine

        Args:
            prog (dict): option to run
            directory (str): directory to run
            common_param (dict): akaikkr inpout parameter
            comment (str): comment to show
            go (str, optional): go keyword. Defaults to "go".
            inputcard (str, optional): inputcard filename. Defaults to None, then inputcard filename is generated automatically.
            outputcard (str, optional): outputcard filename. Defaults to None, then outputcard filename is genarated automatically.
        """
        self.akaikkr_exe = prog["specx"]
        self.args = prog["args"]
        self.no_run = self.args.no_run  # args is not dict, but namespace.
        self.directory = directory
        self.comment = comment
        self.go = go
        self.param = common_param
        self.args = args
        if inputcard is None:
            self.inputcard = _make_inputcard(self.go)
        else:
            self.inputcard = inputcard
        if outputcard is None:
            self.outputcard = _make_outputcard(self.go)
        else:
            self.outputcard = outputcard
        os.makedirs(directory, exist_ok=True)
        msg = " ".join(["###", self._make_comment()])
        print("".join(["-" for char in msg]))
        print(msg)

    @abstractmethod
    def prescript(self, ):
        self.param["go"] = self.go
        self.param["record"] = "init"
        if self.args:
            self.param.update(self.args)

    def get_result_testrun(self, outfile, Awk_k=None):
        """get result for testrun

        Args:
            outfile (str): output filename of akaikkr to analyze
            Awk_k (int, optional): k index to save A(w,k). Defaults to None.

        Returns:
            dict: result summary
        """

        job = AkaikkrJob(self.directory)
        go = job.get_go(outfile)

        dic = {"go": go, "te": job.get_total_energy(outfile),
               "rms": job.get_rms_error(outfile)[-1],
               "tm": job.get_total_moment(outfile),
               "spinlocalmoment": job.get_local_moment(outfile, mode="spin"),
               "orbitallocalmoment": job.get_local_moment(outfile, mode="orbital"),
               "threads": job.get_threads_openmp(outfile),
               "conv": job.check_convergence_go(outfile)}
        if go == "fsm":
            dic.update({"fspin": job.get_fixed_spin_moment(outfile)})
        if go[0] == "j":
            dic.update({"Tc": job.get_curie_temperature(outfile)})
            _jijdf = job.get_jij_as_dataframe(outfile)
            jijselect = _cut_only_nn(_jijdf)
            dic.update({"jij": jijselect})
        if go == "tc":
            dic.update({"Tc": job.get_curie_temperature(outfile)})
        if go == "dos":
            ene, dos_block = job.get_dos_as_list(outfile)
            dos_up = dos_block[0]
            if len(dos_block) > 1:
                dos_dn = dos_block[1]
            else:
                dos_dn = []
            dos_dic = {"energy": ene, "up": dos_up, "dn": dos_dn}
            dic.update({"totalDOS": dos_dic})

            pdos = job.get_pdos_as_dictdataframe(outfile)
            iatom = 0
            l_label = "d"
            updn = "pdos_up"
            df = pdos[updn][iatom]
            ene = df.loc[:, "energy"].values.tolist()
            pdos_up_d = df.loc[:, l_label].values.tolist()
            updn = "pdos_dn"
            pdos_dn_d = []
            if len(pdos[updn]) >= 1:
                df = pdos[updn][iatom]
                pdos_dn_d = df.loc[:, l_label].values.tolist()
            pdos_dic = {"energy": ene, "up": pdos_up_d, "dn": pdos_dn_d,
                        "atom": iatom, "l": l_label}
            dic.update({"pdos": pdos_dic})

        if go == "cnd":
            resistivity = job.get_resistivity(outfile)
            conduct = job.get_conductivity_spin(outfile)
            dic.update({"resis": resistivity, "cnd": conduct})
        if go[:3] == "spc":
            magtyp = job.get_magtyp(outfile)
            pot = job.get_potentialfile(outfile)
            Awk_dic = {}
            if magtyp[0] == "m":  # nspin = 2
                filename = "{}_up.spc".format(pot)
                filepath = os.path.join(job.path_dir, filename)
                Awk = AwkReader(filepath)
                if Awk_k is None:
                    Awk_dic["kpath"] = Awk.kcrt.tolist()
                    Awk_dic["kdist"] = Awk.kdist.tolist()
                    Awk_dic["energy"] = Awk.energy.tolist()
                    Awk_dic["Awk_up"] = Awk.Awk.tolist()
                    filename = "{}_dn.spc".format(pot)
                    filepath = os.path.join(job.path_dir, filename)
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
                    filepath = os.path.join(job.path_dir, filename)
                    Awk = AwkReader(filepath)
                    Awk_dic["Aw_dn"] = Awk.Awk[Awk_k].tolist()
            else:  # nspin = 1
                filename = "{}_up.spc".format(pot)
                filepath = os.path.join(job.path_dir, filename)
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

    def execute(self, displc=False, execute_postscript=True):
        """execute Akaikkr

        Args:
            displc (bool, optional): displc option. Defaults to False. True to add displc
        """
        self.prescript()
        if displc:
            self.param["displc"] = _make_displc(self.param["anclr"])
        job = AkaikkrJob(self.directory)
        inputcard = self.inputcard
        outputcard = self.outputcard
        _param = job.default
        _param.update(self.param)
        if not self.no_run:
            job.make_inputcard(_param, inputcard)
            job.run(self.akaikkr_exe, inputcard, outputcard)

        # total magnetic momenet can be negative when no SCF is requested.
        if self.go == "go":
            totalmoment = job.get_total_moment(outputcard)
            if totalmoment < 0:
                print("  comment: total moment<0. rerun it.")
                _param.update({"record": "2nd"})
                if not self.no_run:
                    job.make_inputcard(_param, inputcard)
                    job.run(self.akaikkr_exe, inputcard, outputcard)

        print("  directory: ", job.path_dir, "  output file: ", outputcard)

        result = self.get_result_testrun(outputcard, Awk_k=0)
        result["comment"] = self._make_comment()
        self.result = result
        if execute_postscript:
            self.postscript()

    @abstractmethod
    def postscript(self, ):
        pass

    def _make_comment(self):
        return self.comment.format(self.go)


class GoFmg(GoGo):
    """run Akaikkr as go="go" after spin flipping
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 flip_list,
                 go="gofmg", **args):
        """[summary]

        Args:

            akaikkr_exe (str): akaikkr program path
            directory (str): directory to run
            common_param (dict): akaikkr inpout parameter
            comment (str): comment to show
            flip_list (list): a list of shortnames to spin flipping. shortnamnes are found in Akaikkr.typeofsites()
            go (str, optional): go keyword. Defaults to "gofmg".
        """
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)

        self.fmg_exe = akaikkr_exe["fmg"]
        self.flip_list = flip_list

    def prescript(self,):
        directory = self.directory
        job = AkaikkrJob(directory)
        outfile = "out_go.log"  # no way to get it
        typeofsites = job.get_type_of_site(outfile)
        potentialfile1, potentialfile2 = job.default["potentialfile"], "pot_fmg.dat"
        if not self.no_run:
            fmg = Fmg(directory)
            fmg.make_inputfile(typeofsites, potentialfile1, potentialfile2,
                               self.flip_list)
            fmg.run(self.fmg_exe)

        self.param["potentialfile"] = potentialfile2
        self.param["go"] = "go"
        self.param["record"] = "2nd"
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self, ):
        pass


class GoDos(GoGo):
    """run Akaikkr as go="dos". Png images are also generated.
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="dos", **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)

    def prescript(self,):
        self.param["go"] = self.go
        self.param["record"] = "2nd"
        self.param["ewidth"] = 2.0  # Default value. It will be updated.
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self):
        plot_dos(self.directory, self.outputcard)
        plot_pdos_all(self.directory, self.outputcard)


class GoTc(GoGo):
    """run Akaikkr as go="tc"
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="tc", **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)

    def prescript(self, ):
        self.param["go"] = self.go
        self.param["record"] = "2nd"
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self,):
        pass


class Goj30(GoGo):
    """run Akaikkr as go="j3.0"
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="j3.0", **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)

    def prescript(self, ):
        self.param["go"] = self.go
        self.param["record"] = "2nd"
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self):
        job = AkaikkrJob(self.directory)
        df = job.get_jij_as_dataframe(self.outputcard)
        filename = "jij.csv"
        filepath = os.path.join(self.directory, filename)
        df.to_csv(filepath, index=False)
        print("  saved to", filepath)
        print()


class Goj(GoGo):
    """run Akaikkr as go="j"
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="j", **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)

    def prescript(self, **args):
        self.param["go"] = self.go
        self.param["record"] = "2nd"
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self):
        job = AkaikkrJob(self.directory)
        df = job.get_jij_as_dataframe(self.outputcard)
        filename = "jij.csv"
        filepath = os.path.join(self.directory, filename)
        df.to_csv(filepath, index=False)
        print("  save to", filepath)
        print()


class GoFsm(GoGo):
    """run Akaikkr as go="fsm"
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="fsm", fspin=1.0, **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)
        self.fspin = fspin

    def prescript(self, **args):
        self.param["go"] = self.go
        self.param["record"] = "init"
        self.param["fspin"] = self.fspin
        self.param["pmix"] = "0.02ch"  # specify Chebyshev mixing
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self):
        pass


class GoCnd(GoGo):
    """run Akaikkr as go="cnd"
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="cnd", **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)

    def prescript(self, **args):
        self.param["go"] = " "+self.go  # extran " " is necessary.
        self.param["record"] = "2nd"
        self.param["ewidth"] = 0.01  # default value.
        self.param["bzqlty"] = 40  # default value. It will be updated later.
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self):
        pass


class GoSpc(GoGo):
    """run Akaikkr as go="spc". Png images are also generated.
    """

    def __init__(self, akaikkr_exe, directory, common_param, comment,
                 go="spc31", go_conv="go", nk=150, **args):
        super().__init__(akaikkr_exe, directory, common_param, comment, go=go, **args)
        self.go_conv = go_conv
        self.nk = nk
        if "fmt" in args:
            self.fmt = args["fmt"]
        else:
            self.fmt = 3
        if "first_connected_kpath" in args:
            self.first_connected_kpath = args["first_connected_kpath"]
        else:
            self.first_connected_kpath = True

    def prescript(self, **args):
        outputcard = "out_{}.log".format(self.go_conv)
        directory = self.directory
        go_job = AkaikkrJob(directory)
        # pymatgen symmetry finder failed in the case of alloys
        # force to change atom name as element.
        struc = go_job.make_pymatgenstructure(
            "out_go.log", change_atom_name=True)

        klabel_filename = KLABEL_FILENAME_
        klabel_filepath = os.path.join(directory, klabel_filename)
        highsymkpath = HighSymKPath(structure=struc,
                                    klabel_filename=klabel_filepath)
        # also make json file if klabel_filename is not None
        kpath = highsymkpath.make_akaikkr_lines(nk=self.nk, fmt=self.fmt,
                                                first_connected_kpath=self.first_connected_kpath)

        _param = {"kpath_raw": kpath}
        self.param.update(_param)
        self.param["go"] = self.go
        self.param["record"] = "2nd"
        args = self.args
        if args:
            self.param.update(args)

    def postscript(self):
        plotband(self.directory)
