# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import sys
import math

import pandas as pd


def _show_df_exist(df):
    def _foo():
        pass

    _df = df.copy()
    for col_name, v in _df.iteritems():
        vlist = []
        for x in v:
            if isinstance(x, type(_foo)):
                vlist.append("-")
            else:
                if math.isnan(x):
                    vlist.append("")
                else:
                    vlist.append("")
        _df[col_name] = vlist
    print(_df)
    print()


def _make_exe_list(df):
    def _foo():
        pass
    exe_list = []
    for index, row in df.iterrows():
        for x in row:
            if isinstance(x, type(_foo)):
                exe_list.append(x)
    return exe_list


class ExeUtil:
    """when a list of function is given, this class shows it as a table.
    """

    def __init__(self, exe_list):
        """initialization routine

        Args:
            exe_list (list): a list of functions
        """
        exe_name_dic = {}
        for exe in exe_list:
            exe_name_dic[exe.__name__] = exe

        # make material name and go name
        materialnames = []
        gonames = []
        for name in exe_name_dic.keys():
            s = name.split("_")
            goname = s[-1]
            materialname = "_".join(s[:-1])
            materialnames.append(materialname)
            gonames.append(goname)

        uniq_materialnames = list(set(materialnames)).sort()
        uniq_gonames = list(set(gonames)).sort()

        df = pd.DataFrame(None, columns=uniq_gonames, index=uniq_materialnames)

        for materialname, goname in zip(materialnames, gonames):
            name = "_".join([materialname, goname])
            df.loc[materialname, goname] = exe_name_dic[name]

        self.df = df

    def exe_list(self):
        """make a list of functions. 
        This shoudle be the same as the exe_list as an input.

        Returns:
            list: a list of functions
        """
        return _make_exe_list(self.df)

    def show_exist(self):
        """show functions as a table.
        """
        _show_df_exist(self.df)
