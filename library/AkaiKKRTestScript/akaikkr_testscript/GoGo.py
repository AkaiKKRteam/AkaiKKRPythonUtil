# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.



import os
from abc import abstractmethod
from pyakaikkr import plotband, KLABEL_FILENAME_, AkaikkrJob, HighSymKPath, plot_dos, plot_pdos_all, Fmg


def _make_displc(anclr, displc=[0, 0, 0]):
    displc_list = []
    for anclr1 in anclr:
        displc1_list = []
        for z in anclr1:
            displc1_list.append(displc)
        displc_list.append(displc1_list)
    return displc_list


def _make_inputcard(go):
    return "inputcard_{}".format(go)


def _make_outputcard(go):
    return "out_{}.log".format(go)


class GoGo:
    """run Akaikkr as go="go"
    """

    def __init__(self, prog,
                 directory,  common_param, comment,
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

        result = job.get_result_testrun(outputcard, Awk_k=0)
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
        self.param["ewidth"] = 2.0
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
        df = job.cut_jij_dataframe(self.outputcard)
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
        df = job.cut_jij_dataframe(self.outputcard)
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
        self.param["ewidth"] = 0.01
        self.param["bzqlty"] = 40
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
