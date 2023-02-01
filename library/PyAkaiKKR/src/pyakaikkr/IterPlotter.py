# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import matplotlib.pyplot as plt
import os
import numpy as np

from .AkaiKkr import AkaikkrJob

from .BasePlotter import BaseEXPlotter


class IterPlotter:
    """plotter for history
    """

    def __init__(self, rms):
        """initialization routine

        Args:
            rms ([float]): history values
        """
        rms = np.array(rms)
        if rms.ndim == 1:
            rms = rms.reshape(1, -1)
        self.rms = rms

    def make(self, output_directory: str, ylabels: list, filename: str,  figsize=(5, 3)):
        """make iteration plot

        Args:
            output_directory (str): output directory
            ylabels (list): ylabels
            filename (str): output filename
            figsize (tuple, optional): figure size. Defaults to (5, 3).
        """
        outputpath = output_directory
        os.makedirs(outputpath, exist_ok=True)
        filepath = os.path.join(outputpath, filename)

        if not isinstance(ylabels, list):
            ylabels = [ylabels]

        fig, axes = plt.subplots(self.rms.shape[0], 1, figsize=figsize)
        if not isinstance(axes, np.ndarray):
            axes = [axes]
        for y, ax, ylabel in zip(self.rms, axes, ylabels):
            x = list(range(len(y)))
            x = np.array(x)
            x += 1
            ax.plot(x, y)
            ax.set_ylabel(ylabel)
            ax.tick_params(axis="x", labelbottom=False)
        # show only the last ticks and labels
        ax.set_xlabel("iteration")
        ax.tick_params(axis="x", labelbottom=True)
        fig.tight_layout()
        fig.savefig(filepath)
        print("saved to", filepath)
        fig.clf()
        plt.close(fig)


class IterEXPlotter(BaseEXPlotter):
    def __init__(self, directory, outfile="out_go.log",  output_directory=None,):
        """
        Args:
            directory (str): directory to save figures
            outfile (str, optional): output filename. Defaults to "out_spc.log".
            pot (str, optional): potential filename. Defaults to "pot.dat".
            output_directory (str, optional): the directory of the output file. Defaults to None.

        """

        super().__init__(directory, outfile, output_directory)

    def make(self, hist_type=["te", "moment", "err"], filename: str = "iter_all.png", figsize=(5, 3)):
        """make history plot from outputfile

        Args:
            hist_type ([str]]): history type te|moment|err. Defauls to ["te", "moment", "err"].
            filename (str): image filename
        """
        job = AkaikkrJob(self.directory)
        rms = []
        for h in hist_type:
            if h == "te":
                value = job.get_te_history(self.outfile)
            elif h == "moment":
                value = job.get_moment_history(self.outfile)
            elif h == "err":
                value = job.get_err_history(self.outfile)
            else:
                raise ValueError("unknown hist_type={}".format(hist_type))
            rms.append(value)

        iterplotter = IterPlotter(rms)
        iterplotter.make(self.output_directory, ylabels=hist_type,
                         filename=filename, figsize=figsize)
