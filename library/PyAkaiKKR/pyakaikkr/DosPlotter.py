# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


import os
import matplotlib.pyplot as plt
import numpy as np
from .AkaiKkr import AkaikkrJob


def plot_dos(directory, outfile, output_direcotry=None,
             yscale="log", figsize=(5, 3)):
    job = AkaikkrJob(directory)
    energy, dos_up, dos_dn = job.cut_dos_float(outfile)
    if len(dos_dn) > 1:  # mag
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


def plot_pdos_all(directory, outfile, output_direcotry=None,
                  yscale="log", figsize=(5, 3)):
    l_label = ["s", "p", "d", "f", "g"]
    job = AkaikkrJob(directory)
    typeofsite = job.get_type_of_site(outfile)

    serial_site = []
    for site in typeofsite:
        for shortname in site["comp_shortname"]:
            serial_site.append(shortname)

    output_format = "spin_separation"
    energy, pdos_block = job.cut_pdos_all(outfile, output_format=output_format)
    energy = np.array(energy)
    # print("energy", energy)
    # print("pdos_block", pdos_block)
    if False:
        # print out format
        if output_format == "spin_separation":
            for pdos_updn in pdos_block:
                for pdos in pdos_updn:
                    pdos = np.array(pdos)
                    print(pdos[:2, :])
        elif output_format == "sequential":
            for pdos in pdos_block:
                pdos = np.array(pdos)
                print(pdos[:2, :])
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


class DosPlotter:
    def __init__(self, directory):
        self.directory = directory

    def show(self, outfile, output_directory=None):
        plot_dos(self.directory, outfile, output_direcotry=output_directory)
        plot_pdos_all(self.directory, outfile,
                      output_direcotry=output_directory)
