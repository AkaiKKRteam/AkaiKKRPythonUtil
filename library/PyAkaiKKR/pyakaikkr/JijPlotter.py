# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import pandas as pd
import matplotlib.pyplot as plt
import os
from .BasePlotter import BaseEXPlotter, BasePlotter
from .AkaiKkr import AkaikkrJob


class JijPlotter(BasePlotter):
    def __init__(self, df, output_directory):
        """initialization routine

        Args:
            df (pd.DataFrame, optional): data. Defaults to None
            output_directory (str): output directory of images.
        """
        super().__init__(output_directory)
        self.df = df

    def make_typepair(self, a=1.0, marker="o",
                      figsize=(5, 3)):
        """plot Jij of all type pairs
        save all combinations of type pairs

        Args:
            a (float, optional): a length. Defaults to 1.0.
            marker (str, optional): matplotlib.plot() maker. Defaults to "o".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """

        df = self.df.copy()
        output_directory = self.output_directory

        xlabel = "distance"
        ylabel = "J_ij(meV)"
        xlabel_fig = "$R$"
        ylabel_fig = "$J_{ij}$ (meV)"

        # plot range
        values = df[xlabel].astype(float).values
        xlim = (values.min(), values.max())
        dx = (xlim[1]-xlim[0])*0.05
        xlim = (xlim[0]-dx, xlim[1]+dx)
        values = df[ylabel].astype(float).values
        ylim = (values.min(), values.max())
        print("debug, ylim", ylim)
        dy = (ylim[1]-ylim[0])*0.05
        ylim = (ylim[0]-dy, ylim[1]+dy)

        # make type pair
        type1 = df["type1"]
        type2 = df["type2"]
        type_pair_list = []
        for t1, t2 in zip(type1, type2):
            s = "-".join([t1, t2])
            type_pair_list.append(s)
        df["typepair"] = type_pair_list
        uniq_type_pair = list(set(type_pair_list))

        # make figures
        for pair in uniq_type_pair:
            _df = df.query("typepair=='{}'".format(pair))

            distance = _df[xlabel].astype(float).values*a
            Jij = _df[ylabel]

            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(distance, Jij, linestyle="-", marker=marker)
            ax.axhline(y=0, linestyle="--", linewidth=1)
            ax.set_xlabel(xlabel_fig)
            ax.set_ylabel(ylabel_fig)
            ax.set_title(pair)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            if True:
                imgfile = "Jij_{}.png".format(pair)
                imgpath = os.path.join(output_directory,
                                       imgfile)
                fig.tight_layout()
                fig.savefig(imgpath)
                print("  saved to", imgpath)
                fig.clf()
                plt.close(fig)

    def make_comppair(self, type1, type2, typeofsite,
                      a=1.0,
                      marker="o",
                      figsize=(5, 3)):
        """plot Jij of specified type1 and type2
        save png images of all the (comp1,comp2) combinations

        Args:
            type1 (str): type1 name
            type2 (str): type2 name
            typeofsite (dict): typeofsite of AkaikkrJob

            a (float, optional): a scale of length. Defaults to 1.0.
            marker (str, optional): matplotlib.plot() marker. Defaults to "o".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        output_directory = self.output_directory
        df = self.df.query("type1=='{}' and type2=='{}'".format(
            type1, type2)).reset_index(drop=True)
        xlabel = "distance"
        ylabel = "J_ij(meV)"
        xlabel_fig = "$R$"
        ylabel_fig = "$J_{ij}$ (meV)"

        # plot range
        values = df[xlabel].astype(float).values*a
        xlim = (values.min(), values.max())
        dx = (xlim[1]-xlim[0])*0.05
        xlim = (xlim[0]-dx, xlim[1]+dx)
        values = df[ylabel].values
        ylim = (values.min(), values.max())
        dy = (ylim[1]-ylim[0])*0.05
        ylim = (ylim[0]-dy, ylim[1]+dy)

        # make comp pair
        comp1 = df["comp1"]
        comp2 = df["comp2"]
        type_pair_list = []
        for t1, t2 in zip(comp1, comp2):
            s = "-".join([str(t1), str(t2)])
            type_pair_list.append(s)
        df_comppair = pd.DataFrame({"comppair": type_pair_list})
        df = pd.concat([df, df_comppair], axis=1)
#        df["comppair"] = type_pair_list # warning occurs
        uniq_type_pair = list(set(type_pair_list))

        shortname_dic = {}
        for component in typeofsite:
            type = component["type"]
            shortname_dic[type] = component["comp_shortname"]

        # make figures

        for pairname in uniq_type_pair:
            s = pairname.split("-")
            comp1, comp2 = int(s[0])-1, int(s[1])-1
            comp1name = shortname_dic[type1][comp1]
            comp2name = shortname_dic[type2][comp2]
            label = "{}-{}".format(comp1name, comp2name)
            _df = df.query("comppair=='{}'".format(pairname))

            distance = _df[xlabel]*a
            Jij = _df[ylabel]

            fig, ax = plt.subplots(figsize=figsize)
            ax.plot(distance, Jij, linestyle="-", marker=marker, label=label)
            ax.axhline(y=0, linewidth=1, linestyle="--")
            ax.set_xlabel(xlabel_fig)
            ax.set_ylabel(ylabel_fig)
            ax.set_title(label)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            if True:
                imgfile = "Jij_{}_{}_{}.png".format(type1, type2, pairname)
                imgpath = os.path.join(output_directory, imgfile)
                fig.tight_layout()
                fig.savefig(imgpath)
                fig.clf()
                plt.close(fig)
                print("  saved to", imgpath)


class JijEXPlotter(BaseEXPlotter, JijPlotter):
    def __init__(self, directory, outfile, output_directory):
        super(BaseEXPlotter, self).__init__(
            directory, outfile, output_directory)

        job = AkaikkrJob(directory)
        df = job.get_jij_as_dataframe()
        super(JijPlotter, self).__init__(df, output_directory)

    # the same members of JijPlotter can be used.
