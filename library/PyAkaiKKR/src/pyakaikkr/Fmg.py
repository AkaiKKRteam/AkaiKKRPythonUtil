# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import os
import subprocess


class Fmg:
    """class to run fmg to flip spin momemts
    """

    def __init__(self, path_dir, ):
        """initialization routine

        Args:
            path_dir (str): directory to run
        """
        self.path_dir = path_dir

    def make_inputfile(self, typeofsites, potentialfile1, potentialfile2,
                       flip_list, input_filename="fmg.input"):
        """make inputfile of fmg

        Args:
            typeofsites (dict): typeofsites output of Akaikkr class
            potentialfile1 (str): input potential file
            potentialfile2 (str): output potential file
            flip_list (list): a list of shortname, which can be found in typeofsites
            input_filename (str, optional): input filename of fmg. Defaults to "fmg.input".
        """
        self.input_filename = input_filename

        # generate serial list
        comp_shortname_list = []
        for component in typeofsites:
            for comp_shortname in component["comp_shortname"]:
                comp_shortname_list.append(comp_shortname)

        # initial values
        all_flip_dic = {}
        for comp in comp_shortname_list:
            all_flip_dic[comp] = 0

        # flip if it exists
        for comp in flip_list:
            all_flip_dic[comp] = 1

        ilist1 = [i+1 for i in range(len(comp_shortname_list))]
        ilist1 = list(map(str, ilist1))
        ilist2 = []
        for icmp, comp in enumerate(comp_shortname_list):
            i = all_flip_dic[comp]
            ilist2.append((icmp+1)*(-1)**i)
        ilist2 = list(map(str, ilist2))

        # generate an input file
        lines = []
        lines.append("{} {}".format(potentialfile1, " ".join(ilist1)))
        lines.append("{} {}".format(potentialfile2, " ".join(ilist2)))
        fmg_inputpath = os.path.join(self.path_dir, self.input_filename)
        with open(fmg_inputpath, "w") as f:
            f.write("\n".join(lines))
        print("  save to", fmg_inputpath)

    def run(self, fmg_exe, outfile="fmg.output"):
        """run fmg

        Args:
            fmg_exe (str): a path to fmg
            outfile (str, optional): output filename of fmg. Defaults to "fmg.output".
        """
        cmd = "cd {}; {} < {} > {}".format(
            self.path_dir, fmg_exe,  self.input_filename, outfile)
        ret = subprocess.call(cmd, shell=True)
        return ret
