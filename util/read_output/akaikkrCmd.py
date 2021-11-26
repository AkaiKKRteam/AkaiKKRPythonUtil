# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from collections import OrderedDict
import os
import sys
from typing import OrderedDict
import click

__verson__ = '2021.10.18'

_path = os.path.split(os.path.split(__file__)[0])[0]
sys.path.append(_path)
print("temporary: add python path", _path)

if "pyakaikkr" not in sys.modules:
    from pyakaikkr import *

_OUTPUTPATH_HELP_ = "output directory"
_PNG_HELP_ = "create png file"
_FILE_HELP_ = "create output file"
_OUTPUTTYPE_HELP_ = "choose output type"


def serial_shortname(tos):
    shortnames = []
    for site in tos:
        shortnames.extend(site["comp_shortname"])
    return shortnames


@click.group()
def cmd():
    """read AkaiKKR output file
    """
    pass


@cmd.command()
def converged():
    click.echo('converged')


@cmd.command()
@click.argument("filename", type=str, required=True)
def converged(filename):
    """Extract converged or not"""
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    v = job.check_convergence_go(filename)
    print(sys._getframe().f_code.co_name, v)


@cmd.command()
@click.argument("filename", type=str, required=True)
def lasterr(filename):
    """Extract the last err"""
    job = AkaikkrJob(".")
    v = job.get_rms_error(filename)
    print(sys._getframe().f_code.co_name, v[-1])


@cmd.command()
@click.argument("filename", type=str, required=True)
def totalenergy(filename):
    """Extract total energy"""
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    v = job.get_total_energy(filename)
    print(sys._getframe().f_code.co_name, v)


@cmd.command()
@click.argument("filename", type=str, required=True)
def totalmoment(filename):
    """Extract total moment"""
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    v = job.get_total_moment(filename)
    print(sys._getframe().f_code.co_name, v)


def _spinlocalmoment(filename):
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    shortnames = serial_shortname(tos)
    tm = job.get_local_moment(filename)
    print(sys._getframe().f_code.co_name,)
    d = []
    for name, v in zip(shortnames, tm):
        d.append((name, v))
    print(d)


@cmd.command()
@click.argument("filename", type=str, required=True)
def localmoment(filename):
    """Extract spin local moment, the same as spinlocalmoment"""
    _spinlocalmoment(filename)


@cmd.command()
@click.argument("filename", type=str, required=True)
def spinlocalmoment(filename):
    """Extract spin local moment"""
    _spinlocalmoment(filename)


@cmd.command()
@click.argument("filename", type=str, required=True)
def orbitallocalmoment(filename):
    """Extract orbital local moment"""
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    shortnames = serial_shortname(tos)
    tm = job.get_local_moment(filename, mode="orbital")
    print(sys._getframe().f_code.co_name,)
    d = []
    for name, v in zip(shortnames, tm):
        d.append((name, v))
    print(d)


@cmd.command()
@click.argument("filename", type=str, required=True)
def charge(filename):
    """Extract charge"""
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    shortnames = serial_shortname(tos)
    tm = job.get_type_charge(filename)
    print(sys._getframe().f_code.co_name,)
    d = []
    for name, v in zip(shortnames, tm):
        d.append((name, v))
    print(d)


@cmd.group()
def hist():
    """Extract te/moment/err histories"""
    pass


@hist.command("all")
@click.argument("filename", type=str, required=True)
@click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
              default='output', show_default=True)
@click.option('--outputtype', '-t', help=_OUTPUTTYPE_HELP_,
              type=click.Choice(["png", "csv", "console"]), default="png")
def histall(filename, outputpath, outputtype):
    """Extract all histories"""
    job = AkaikkrJob(".")
    m = job.get_moment_history(filename)
    te = job.get_te_history(filename)
    err = job.get_err_history(filename)
    xs = [m, te, err]
    labels = ["moment", "te", "err"]
    if outputtype == "csv":
        df = pd.DataFrame({"moment": m, "te": te, "err": err})
        os.makedirs(outputpath, exist_ok=True)
        filepath = os.path.join(outputpath, "hist.csv")
        df.to_csv(filepath, index=False)
        print("saved to", filepath)
    elif outputtype == "console":
        for v, label in zip([m, te, err], labels):
            print(label, v)
    elif outputtype == "png":
        iterplotter = IterPlotter(xs)
        os.makedirs(outputpath, exist_ok=True)
        iterplotter.show(outputpath, filename="hist.png",
                         ylabels=labels, figsize=(5, 8))


@ hist.command("moment")
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
               default='output', show_default=True)
@click.option('--outputtype', '-t', help=_OUTPUTTYPE_HELP_,
              type=click.Choice(["png", "csv", "console"]), default="png")
def momenthist(filename, outputpath, outputtype):
    """Extract moment history"""
    job = AkaikkrJob(".")
    tm = job.get_moment_history(filename)
    if outputtype == "console":
        print(sys._getframe().f_code.co_name, tm)
    elif outputtype == "csv":
        df = pd.DataFrame({"tm": tm})
        os.makedirs(outputpath, exist_ok=True)
        filepath = os.path.join(outputpath, "momenthist.csv")
        df.to_csv(filepath, index=False)
        print("saved to", filepath)
    elif outputtype == "png":
        os.makedirs(outputpath, exist_ok=True)
        rmsplotter = IterPlotter(tm)
        rmsplotter.show(outputpath, "moment", "momenthist.png")


@ hist.command("err")
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
               default='output', show_default=True)
@click.option('--outputtype', '-t', help=_OUTPUTTYPE_HELP_,
              type=click.Choice(["png", "csv", "console"]), default="png")
def errhist(filename, outputpath, outputtype):
    """Extract err history"""
    job = AkaikkrJob(".")
    tm = job.get_err_history(filename)
    if outputtype == "console":
        print(sys._getframe().f_code.co_name, tm)
    elif outputtype == "csv":
        df = pd.DataFrame({"tm": tm})
        os.makedirs(outputpath, exist_ok=True)
        filepath = os.path.join(outputpath, "errhist.csv")
        df.to_csv(filepath, index=False)
        print("saved to", filepath)
    elif outputtype == "png":
        os.makedirs(outputpath, exist_ok=True)
        rmsplotter = IterPlotter(tm)
        rmsplotter.show(outputpath, "err", "errhist.png")


@ hist.command("te")
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
               default='output', show_default=True)
@click.option('--outputtype', '-t', help=_OUTPUTTYPE_HELP_,
              type=click.Choice(["png", "csv", "console"]), default="png")
def tehist(filename, outputpath, outputtype):
    """Extract te history"""
    job = AkaikkrJob(".")
    tm = job.get_te_history(filename)
    if outputtype == "console":
        print(sys._getframe().f_code.co_name, tm)
    elif outputtype == "csv":
        df = pd.DataFrame({"tm": tm})
        os.makedirs(outputpath, exist_ok=True)
        filepath = os.path.join(outputpath, "tehist.csv")
        df.to_csv(filepath, index=False)
        print("saved to", filepath)
    elif outputtype == "png":
        os.makedirs(outputpath, exist_ok=True)
        rmsplotter = IterPlotter(tm)
        rmsplotter.show(outputpath, "te", "tehist.png")


@ cmd.command()
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
               default='output', show_default=True)
@ click.option("--format", "-f", help="format type",
               type=click.Choice(["xsf", "cif", "poscar"]),
               default="cif")
@ click.option('--outputtype', '-t', help=_OUTPUTTYPE_HELP_,
               type=click.Choice(["file", "console"]),
               default="file")
def structure(filename, outputpath, format, outputtype):
    """Extract lattice vector and atomic positions"""
    job = AkaikkrJob(".")
    struc = job.make_pymatgenstructure(filename, unit="ang")
    writer = StructureWriter(struc)
    os.makedirs(outputpath, exist_ok=True)
    if format == "xsf":
        outputfile = os.path.join(outputpath, "structure.xsf")
        if outputtype == "console":
            outputfile = None
        writer.xsf(outputfilename=outputfile)
    elif format == "cif":
        outputfile = os.path.join(outputpath, "structure.cif")
        if outputtype == "console":
            outputfile = None
        writer.cif(outputfilename=outputfile)
    elif format == "poscar":
        outputfile = os.path.join(outputpath, "POSCAR.vasp")
        if outputtype == "console":
            outputfile = None
        writer.poscar(outputfilename=outputfile)
    else:
        print("unknown format", outputtype)
        raise ValueError("unknown format")


@ cmd.command()
@ click.argument("filename", type=str, required=True)
def typeofsite(filename,):
    """Extract type of site"""
    job = AkaikkrJob(".")
    tos = job.get_type_of_site(filename)
    print(sys._getframe().f_code.co_name,)
    for t in tos:
        print(t)


@ cmd.command()
@ click.argument("filename", type=str, required=True)
def edelt(filename,):
    """Extract edelt"""
    job = AkaikkrJob(".")
    v = job.get_edelt(filename)
    print(sys._getframe().f_code.co_name, v)


@ cmd.command()
@ click.argument("filename", type=str, required=True)
def ewidth(filename,):
    """Extract ewidth
    """
    job = AkaikkrJob(".")
    v = job.get_ewidth(filename)
    print(sys._getframe().f_code.co_name, v)


@ cmd.command()
@ click.argument("filename", type=str, required=True)
def fspin(filename,):
    """Extract fixed spin moment
    """
    job = AkaikkrJob(".")
    v = job.get_fixed_spin_moment(filename)
    print(sys._getframe().f_code.co_name, v)


@ cmd.command()
@ click.argument("filename", type=str, required=True)
def Tc(filename,):
    """Extract Curie temperature
    """
    job = AkaikkrJob(".")
    v = job.get_resistivity(filename)
    print(sys._getframe().f_code.co_name, v)


@ cmd.command("R")
@ click.argument("filename", type=str, required=True)
def R(filename,):
    """Extract Resistivity and Conductivity per spin
    """
    job = AkaikkrJob(".")
    v = job.get_resistivity(filename)
    print("R", v)
    v = job.get_conductivity_spin(filename)
    print("cnd", v)


@ cmd.command("Jij")
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o',
               help=_OUTPUTPATH_HELP_, default='output', show_default=True)
@ click.option('--outputtype', '-t', help=_OUTPUTTYPE_HELP_,
               type=click.Choice(["png", "csv", "console"]), default="png")
def Jij(filename, outputpath, outputtype):
    """Generate Jij plot/csv
    """
    job = AkaikkrJob(".")
    os.makedirs(outputpath, exist_ok=True)

    if outputtype == "csv":
        df = job.cut_jij_dataframe(filename)
        jij_filename = os.path.join(outputpath, "Jij.csv")
        df.to_csv(jij_filename, index=False)
        print(jij_filename, "is made.")
    elif outputtype == "console":
        df = job.cut_jij_dataframe(filename)
        print(df)
    elif outputtype == "png":
        typeofsite = job.get_type_of_site(filename)
        df = job.cut_jij_dataframe(filename)
        jij_filename = os.path.join(outputpath, "Jij.csv")
        df.to_csv(jij_filename, index=False)
        print(jij_filename, "is made.")
        jijplotter = JijPlotter(outputpath, filename="Jij.csv")
        jijplotter.plot_typepair(output_directory=outputpath)
        typepairs = []
        for type1, type2 in zip(df["type1"], df["type2"]):
            typepairs.append("{}-{}".format(type1, type2))
        typepairs = list(set(typepairs))
        for typepair in typepairs:
            types = typepair.split("-")
            jijplotter.plot_comppair(
                types[0], types[1], typeofsite,
                output_directory=outputpath)


@ cmd.command()
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
               default='output', show_default=True)
def dos(filename, outputpath,):
    """Generate DOS and PDOS plot
    """
    job = AkaikkrJob(".")
    os.makedirs(outputpath, exist_ok=True)
    dosplot = DosPlotter(".")
    dosplot.show(filename, outputpath)


@ cmd.command()
@ click.argument("filename", type=str, required=True)
@ click.option('--outputpath', '-o', help=_OUTPUTPATH_HELP_,
               default='output', show_default=True)
@ click.option('--klabelfile',
               help="klabelfile must in the same directory as FILENMAE",
               default=None)
def spc(filename, outputpath, klabelfile):
    """Generate A(w,k) plot
    """
    s = os.path.split(filename)
    dir_part = os.path.join(*s[:-1])
    filename_ = s[-1]
    print(dir_part, filename_)

    plotband(dir_part, outfile=filename_,
             output_directory=outputpath,
             klabel_filename=klabelfile)


if __name__ == '__main__':
    cmd()
