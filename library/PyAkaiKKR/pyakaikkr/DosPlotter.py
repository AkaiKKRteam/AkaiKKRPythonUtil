# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


from .BasePlotter import BasePlotter
import os
import matplotlib.pyplot as plt
import numpy as np

from .BasePlotter import BaseEXPlotter
from .AkaiKkr import AkaikkrJob


def _plot_dos(energy, dos_block, output_direcotry=None,
              yscale="log", figsize=(5, 3)):
    if len(dos_block) > 1:  # mag
        dos_up, dos_dn = dos_block[0], dos_block[1]
        _figsize = (figsize[0]*2, figsize[1])
        fig, axes = plt.subplots(1, 2, figsize=_figsize)
        ax = axes[0]
        ax.plot(energy, dos_up)
        ax.set_title("up")
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")
        ax = axes[1]
        ax.plot(energy, dos_dn)
        ax.set_title("dn")
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")

    else:  # monmag
        dos_up = dos_block[0]
        _figsize = figsize
        fig, ax = plt.subplots(figsize=_figsize)
        ax.plot(energy, dos_up)
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")

    if output_direcotry is None:
        output_direcotry = directory
    fig.tight_layout()
    imgfile = "dos.png"
    imgfilepath = os.path.join(output_direcotry, imgfile)
    fig.savefig(imgfilepath)
    print("  saved to", imgfilepath)
    # fig.show()
    fig.clf()
    plt.close(fig)


def _EX_plot_dos(directory, outfile, output_direcotry=None,
                 yscale="log", figsize=(5, 3)):
    job = AkaikkrJob(directory)
    energy, dos_block = job.get_dos_as_list(outfile)

    if len(dos_block) > 1:  # mag
        dos_up, dos_dn = dos_block[0], dos_block[1]
        _figsize = (figsize[0]*2, figsize[1])
        fig, axes = plt.subplots(1, 2, figsize=_figsize)
        ax = axes[0]
        ax.plot(energy, dos_up)
        ax.set_title("up")
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")
        ax = axes[1]
        ax.plot(energy, dos_dn)
        ax.set_title("dn")
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")

    else:  # monmag
        dos_up = dos_block[0]
        _figsize = figsize
        fig, ax = plt.subplots(figsize=_figsize)
        ax.plot(energy, dos_up)
        ax.set_yscale(yscale)
        ax.set_xlabel("E(Ry)-EF")
        ax.set_ylabel("DOS")

    if output_direcotry is None:
        output_direcotry = directory
    fig.tight_layout()
    imgfile = "dos.png"
    imgfilepath = os.path.join(output_direcotry, imgfile)
    fig.savefig(imgfilepath)
    print("  saved to", imgfilepath)
    # fig.show()
    fig.clf()
    plt.close(fig)


def _plot_pdos_all(energy, pdos_block, typeofsite, output_direcotry=None,
                   yscale="log", figsize=(5, 3)):
    l_label = ["s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m"]

    serial_site = []
    for site in typeofsite:
        for shortname in site["comp_shortname"]:
            serial_site.append(shortname)

    energy = np.array(energy)

    if output_direcotry is None:
        output_direcotry = directory

    if len(pdos_block) > 1:  # up and down
        pdos_up, pdos_dn = pdos_block[0], pdos_block[1]

        _figsize = (figsize[0]*2, figsize[1])
        for icmp, (up_atom, dn_atom, title) in enumerate(zip(pdos_up, pdos_dn, serial_site)):
            fig, axes = plt.subplots(1, 2, figsize=_figsize)

            up_atom = np.array(up_atom)
            dn_atom = np.array(dn_atom)
            ax = axes[0]
            for l in range(up_atom.shape[1]):
                ax.plot(energy, up_atom[:, l], label=l_label[l])
            ax.legend()
            ax.set_title("up")
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            ax = axes[1]
            for l in range(dn_atom.shape[1]):
                ax.plot(energy, dn_atom[:, l], label=l_label[l])
            ax.legend()
            ax.set_title("dn")
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            fig.suptitle(title)
            fig.tight_layout()
            imgfile = "pdos_{}.png".format(icmp)
            imgfilepath = os.path.join(output_direcotry, imgfile)
            fig.savefig(imgfilepath)
            print("  saved to", imgfilepath)
            # fig.show()
            fig.clf()
            plt.close(fig)
        print()

    else:  # up only
        pdos_up = pdos_block[0]

        _figsize = (figsize[0], figsize[1])
        for icmp, (up_atom, title) in enumerate(zip(pdos_up, serial_site)):
            fig, ax = plt.subplots(1, 1, figsize=_figsize)

            up_atom = np.array(up_atom)
            for l in range(up_atom.shape[1]):
                ax.plot(energy, up_atom[:, l], label=l_label[l])
            ax.legend()
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            ax.set_title(title)
            fig.tight_layout()
            imgfile = "pdos_{}.png".format(icmp)
            imgfilepath = os.path.join(output_direcotry, imgfile)
            fig.savefig(imgfilepath)
            print("  saved to", imgfilepath)
            # fig.show()
            fig.clf()
            plt.close(fig)
        print()


def _EX_plot_pdos_all(directory, outfile, output_direcotry=None,
                      yscale="log", figsize=(5, 3)):
    l_label = ["s", "p", "d", "f", "g", "h", "i", "j", "k", "l", "m"]
    job = AkaikkrJob(directory)
    typeofsite = job.get_type_of_site(outfile)

    serial_site = []
    for site in typeofsite:
        for shortname in site["comp_shortname"]:
            serial_site.append(shortname)

    output_format = "spin_separation"
    energy, pdos_block = job.get_pdos_as_list(
        outfile, output_format=output_format)
    energy = np.array(energy)

    if output_direcotry is None:
        output_direcotry = directory

    if len(pdos_block) > 1:  # up and down
        pdos_up, pdos_dn = pdos_block[0], pdos_block[1]

        _figsize = (figsize[0]*2, figsize[1])
        for icmp, (up_atom, dn_atom, title) in enumerate(zip(pdos_up, pdos_dn, serial_site)):
            fig, axes = plt.subplots(1, 2, figsize=_figsize)

            up_atom = np.array(up_atom)
            dn_atom = np.array(dn_atom)
            ax = axes[0]
            for l in range(up_atom.shape[1]):
                ax.plot(energy, up_atom[:, l], label=l_label[l])
            ax.legend()
            ax.set_title("up")
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            ax = axes[1]
            for l in range(dn_atom.shape[1]):
                ax.plot(energy, dn_atom[:, l], label=l_label[l])
            ax.legend()
            ax.set_title("dn")
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            fig.suptitle(title)
            fig.tight_layout()
            imgfile = "pdos_{}.png".format(icmp)
            imgfilepath = os.path.join(output_direcotry, imgfile)
            fig.savefig(imgfilepath)
            print("  saved to", imgfilepath)
            # fig.show()
            fig.clf()
            plt.close(fig)
        print()

    else:  # up only
        pdos_up = pdos_block[0]

        _figsize = (figsize[0], figsize[1])
        for icmp, (up_atom, title) in enumerate(zip(pdos_up, serial_site)):
            fig, ax = plt.subplots(1, 1, figsize=_figsize)

            up_atom = np.array(up_atom)
            for l in range(up_atom.shape[1]):
                ax.plot(energy, up_atom[:, l], label=l_label[l])
            ax.legend()
            ax.set_yscale(yscale)
            ax.set_xlabel("E(Ry)-EF")
            ax.set_ylabel("PDOS")
            ax.set_title(title)
            fig.tight_layout()
            imgfile = "pdos_{}.png".format(icmp)
            imgfilepath = os.path.join(output_direcotry, imgfile)
            fig.savefig(imgfilepath)
            print("  saved to", imgfilepath)
            # fig.show()
            fig.clf()
            plt.close(fig)
        print()


class PDosPlotter(BasePlotter):
    def __init__(output_directory):
        super.__init__(output_directory)

    def make(self, energy, pdos_block, typeofsites,
             yscale="log", figsize=(5, 3)):
        """make PDOS from output of akaikkr.get_pdos_as_list().

        energy, pdos_block is output of  akaikkr.get_pdos_as_list().
        typeofsites is output of akaikkr.get_type_of_site().

        Args:
            energy ([float]): energy.
            pdos_block ([[float]]): pdos_block.
            typeofsites (list): type of sites.
            yscale (str, optional): y scale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        output_directory = self.output_directory
        _plot_pdos_all(energy, pdos_block, typeofsites,
                       output_direcotry=self.output_directory,
                       yscale=yscale, figsize=figsize)


class PDosEXPlotter(BaseEXPlotter):
    def __init__(self, directory, outfile, output_directory):
        """
        Args:
            directory (str): directory to save figures
            outfile (str, optional): output filename. Defaults to "out_go.log".
            output_directory (str, optional): output directory Defaults to None.
        """
        super().__init__(directory, outfile, output_directory)

    def make(self, output_directory=None, yscale="log", figsize=(5, 3)):
        """make PDOS from output file.

        Args:
            output_directory (str, optional): output directory of the image file. Defaults to None.
            yscale (str, optional): y scale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        _EX_plot_pdos_all(self.directory, self.outfile,
                          output_direcotry=output_directory, yscale=yscale,
                          figsize=figsize)


class DosPlotter(BasePlotter):
    def __init__(self, output_directory):
        super().__init__(output_directory)

    def make(self, energy, dos_block, yscale="log",
             figsize=(5, 3)):
        """make DOS from akaikkr.get_dos_as_list()

        enegy and dos_block is output of akaikkr.get\dos_as_list()
        Args:
            energy ([float]): energy mesh.
            dos_block ([[float]]): dos block.
            yscale (str, optional): y scale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        _plot_dos(energy, dos_block,
                  output_direcotry=self.output_directory, yscale=yscale,
                  figsize=figsize)


class DosEXPlotter(BaseEXPlotter):
    def __init__(self, directory, outfile, output_directory):
        """
        Args:
            directory (str): directory to save figures
            outfile (str): output filename. 
            output_directory (str): output directory
        """
        super().__init__(directory, outfile, output_directory)

    def make(self, output_directory=None, yscale="log", figsize=(5, 3)):
        """make DOS from outputfile

        Args:
            output_directory (str, optional): output directory of the image file. Defaults to None.
            yscale (str, optional): yscale. Defaults to "log".
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        _EX_plot_dos(self.directory, self.outfile,
                     output_direcotry=self.output_directory, yscale=yscale,
                     figsize=figsize)
