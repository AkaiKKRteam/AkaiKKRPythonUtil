# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import numpy as np
from collections import OrderedDict
import pandas as pd


class ResultUtil:
    """from a dict of result. This class shows it as a table.
    """

    def __init__(self, result_dic):
        """initialization routine.

        Args:
            result_dic (dict): result as a dict
        """

        check_dic = result_dic
        materialnames = []
        gonames = []
        for name in check_dic.keys():
            s = name.split("_")
            goname = s[-1]
            materialname = "_".join(s[:-1])
            materialnames.append(materialname)
            gonames.append(goname)

        uniq_materialnames = list(set(materialnames))
        uniq_materialnames.sort()
        uniq_gonames = list(set(gonames))
        uniq_gonames.sort()

        content = []
        for row in uniq_materialnames:
            v = [np.nan for i in uniq_gonames]
            content.append(v)

        df = pd.DataFrame(content, columns=uniq_gonames,
                          index=uniq_materialnames)
        for materialname, goname in zip(materialnames, gonames):
            name = "_".join([materialname, goname])
            df.loc[materialname, goname] = check_dic[name]
        df = df.fillna("")
        self.df = df

    def show(self):
        print(self.df)
