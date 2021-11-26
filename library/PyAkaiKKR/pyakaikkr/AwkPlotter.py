# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from pymatgen.core.periodic_table import Element
import matplotlib.pyplot as plt
import os
import json
from .AkaiKkr import AkaikkrJob
from .AwkReader import AwkReader

KLABEL_FILENAME_ = "klabel.json"


class AwkPlotter(AwkReader):
    """Plot A(w,k)
    """

    def __init__(self, filepath):
        """initialization routine

        Args:
            filepath (str): filepath of Akaikkr output filename run as go="spc*"
        """

        super().__init__(filepath)
        if False:
            self.kcrt, self.Awk, self.kdist, self.energy = _load_spc_format_type_a(
                filepath)

    def show(self, klabel=None, ylabel="$E(Ry)-E_F$",
             ax=None, figsize=None, imgfilename=None,
             show_kgrid=True, show_ef=True, show_colorbar=False,
             linecolor="w", linewidth=2,
             **pcolorargs):
        """show A(w,k)

        also make imagefile if imgfilename is not None

        Args:
            klabel ([str], optional): a list of klabel. Defaults to None.
            ylabel (str, optional): ylabel string. Defaults to "(Ry)-E_F$".
            ax (matplotlib.axes, optional): axes of matplotlib. Defaults to None.
            figsize ((int,int), optional): figsize. Defaults to None.
            imgfilename (str, optional): image filename. Defaults to None.
            show_kgrid (bool, optional): show k grid line. Defaults to True.
            show_ef (bool, optional): show ef line. Defaults to True.
            show_colorbar (bool, optional): show colorbar. Defaults to False.
            linecolor (str, optional): line color of grids and ef. Defaults to "w".
            linewidth (int, optional): line width of grid and ef. Defaults to 2.
        """             
        kdist = self.kdist
        Awk = self.Awk
        energy = self.energy
        kcrt = self.kcrt
        fig = None
        if ax is None:
            fig, ax = plt.subplots(figsize=figsize)
            fig.tight_layout()

        plt.rcParams['pcolor.shading'] = 'auto'
        self.pcolormesh = ax.pcolormesh(
            kdist, energy, Awk.T, **pcolorargs,)
        ax.set_ylabel(ylabel)
        ax.set_xticks(kdist[kcrt])

        if klabel is not None:
            ax.set_xticklabels(klabel)
        else:
            ax.get_xaxis().set_visible(False)
        if show_kgrid:
            for x in kdist[kcrt]:
                ax.axvline(x, c=linecolor, lw=linewidth)
        if show_ef:
            ax.axhline(0.0, c=linecolor, lw=linewidth)

        if fig and show_colorbar:
            colorbar = fig.colorbar(
                self.pcolormesh, ax=ax, orientation="vertical")

        if fig and imgfilename is not None:
            fig.savefig(imgfilename)

        if fig:
            fig.clf()
            plt.close(fig)


def plotband(directory, pot="pot.dat", outfile="out_go.log",
             output_directory=None, klabel_filename=KLABEL_FILENAME_):
    """plot band A(w,k).
    band files are {potetialfile}_up.spc and {potetialfile}_dn.spc

    Args:
        directory (str): directory to save figures
        pot (str, optional): potential filename. Defaults to "pot.dat".
        outfile (str, optional): output filename. Defaults to "out_go.log".
        output_directory (str, optional): the directory of the output file. Defaults to None.
        klabel_filename (str, optional): klabel filename. Defaults to KLABEL_FILENAME_.
    """
    if output_directory is None:
        output_directory = directory
    job = AkaikkrJob(directory)
    magtyp = job.get_magtyp(outfile)
    param_default = job.default
    pot = job.get_potentialfile(outfile)

    klabel = None
    if klabel_filename is not None:
        kpath_filepath = os.path.join(directory, klabel_filename)
        with open(kpath_filepath) as f:
            content = json.load(f)
            print("debug: file, loaded content",kpath_filepath, content)
            klabel = []
            for vpath in content["kpath"]:
                for v in vpath:
                    k = str(list(v.keys())[0])
                    k = "$"+k+r"$"
                    klabel.append(k)

    filename = "{}_up.spc".format(pot)
    filepath = os.path.join(directory, filename)
    Awk = AwkPlotter(filepath)
    pcolorargs = {"cmap": "magma"}
    imgfilename = os.path.join(output_directory, "Awk_up.png")
    os.makedirs(output_directory, exist_ok=True)
    Awk.show(klabel=klabel, linecolor="green", show_colorbar=True,
             imgfilename=imgfilename, **pcolorargs)
    print("  saved to", imgfilename)

    if magtyp[0] == "m":  # mag
        filename = "{}_dn.spc".format(pot)
        filepath = os.path.join(directory, filename)
        Awk = AwkPlotter(filepath)
        pcolorargs = {"cmap": "magma"}
        imgfilename = os.path.join(output_directory, "Awk_dn.png")
        Awk.show(klabel=klabel, linecolor="green", show_colorbar=True,
                 imgfilename=imgfilename, **pcolorargs)
        print("  saved to", imgfilename)

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        filenames = ["{}_up.spc".format(pot), "{}_dn.spc".format(pot)]
        for filename, ax in zip(filenames, axes):
            filepath = os.path.join(directory, filename)
            Awk = AwkPlotter(filepath)
            pcolorargs = {"cmap": "magma", "alpha": 0.9}
            tickparams = {"labelsize": 12, "labelcolor": "blue"}
            Awk.show(klabel, ax=ax, linewidth=0, **pcolorargs)
            ax.tick_params("both", **tickparams)
            ax.set_ylabel("$E$(Ry)", c="blue")
        imgfilename = os.path.join(output_directory, "Awk_both.png")
        fig.savefig(imgfilename)
        fig.clf()
        plt.close(fig)
        print("  saved to", imgfilename)
    print()
