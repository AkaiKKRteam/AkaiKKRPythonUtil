# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import matplotlib.pyplot as plt
import os
import numpy as np


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
