# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import os


class BasePlotter:
    def __init__(self, output_directory=".", create_directory=True):
        self.output_directory = output_directory


class BaseEXPlotter:
    def __init__(self, directory, outfile, output_directory=None, create_directory=True):
        self.directory = directory
        self.outfile = outfile
        if output_directory is None:
            self.output_directory = directory
        else:
            self.output_directory = output_directory
        if create_directory:
            os.makedirs(output_directory, exist_ok=True)
