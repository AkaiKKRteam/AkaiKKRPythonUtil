# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import json
import os
from collections import OrderedDict
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error
import click


def make_filepathlist_two(directory1, directory2, kkrtype, resultjson):

    resultjson = resultjson
    directorys = [directory1, directory2]
    kkrtype = kkrtype
    idlist = []
    filepathlist = []
    for i, dir_ in enumerate(directorys):
        filepath = os.path.join(dir_, kkrtype, "tests", resultjson)
        filepathlist.append(filepath)
        idlist.append(str(i))
    print(filepathlist)
    print(idlist)
    return filepathlist, idlist


def make_filepathlist_history(prefix, datelist, compiler, kkrtype, resultjson):
    compiler = compiler
    prefix = prefix
    akaikkrtype = kkrtype
    datelist = datelist
    resultjson = resultjson

    filepathlist = []
    idlist = []
    for date in datelist:
        parentdir = ".".join([prefix, date, compiler])
        filepath = os.path.join(parentdir, akaikkrtype, "tests", resultjson)
        filepathlist.append(filepath)
        idlist.append(date)
    print(filepathlist)
    print(idlist)
    return filepathlist, idlist


def compare_them(go, target, subtarget, **ctx):
    metric = ctx["metric"]
    show_df = ctx["show_df"]
    filepathlist = ctx["filepathlist"]
    idlist = ctx["idlist"]
    results = OrderedDict()
    allvaluesdic = {}
    dateprocesslist = []

    for filepath, date in zip(filepathlist, idlist):
        print("read",filepath)
        with open(filepath) as f:
            dic = json.load(f)
            dic = dic["result"]

            results[date] = dic
            valueslist = []
            for materialname, materialcontent in dic.items():
                if go == "jij":
                    # search j3.0 or j
                    flag = False
                    if materialname.endswith("_"+"j3.0"):
                        originalmaterialname = materialname
                        materialname = materialname.replace("_j3.0", "_jij")
                        flag = True
                    elif materialname.endswith("_"+"j"):
                        originalmaterialname = materialname
                        materialname = materialname.replace("_j", "_jij")
                        flag = True
                else:
                    flag = materialname.endswith("_"+go)
                if flag:
                    if target == "totalDOS" or target == "pdos":
                        if subtarget in materialcontent[target]:
                            value = materialcontent[target][subtarget]
                        else:
                            value = None
                    elif target == "jij":
                        value = materialcontent[target]
                        df = pd.DataFrame(value[1:], columns=value[0])
                        value = df["J_ij(meV)"].astype(float).values
                    elif target == "Awk":
                        if subtarget in materialcontent[target]:
                            value = materialcontent[target][subtarget]
                        else:
                            value = None
                    else:
                        value = materialcontent[target]

                    if value is not None:
                        if (isinstance(value, list) and len(value) >= 1) or \
                           (isinstance(value, np.ndarray) and value.shape[0] >= 1) or \
                                isinstance(value, float):
                            if materialname not in allvaluesdic:
                                # print("initialize list", materialname)
                                allvaluesdic[materialname] = []
                            # value = np.array(value)
                            allvaluesdic[materialname].append(value)
                            dateprocesslist.append(date)
                        else:
                            print("no value in", materialname,
                                  target, subtarget, )
                    else:
                        print("no value in", materialname, target, subtarget)

    dateprocesslist = list(set(dateprocesslist))
    dateprocesslist.sort()
    df = pd.DataFrame(allvaluesdic, index=dateprocesslist)

    if show_df:
        print(df.T)
    print("index", df.index.tolist())
    print("columns", df.columns.tolist())

    changed = False
    for col in df.columns:
        vall = df[col].values
        flag = True
        for v in vall[1:]:
            if not np.all(vall[0] == v):
                flag = False
        if not flag:
            print("values changed", col)
            try:
                v = float(vall[0])
                type_ = "scalar"
            except TypeError:
                type_ = "array"
            if type_ == "scalar":
                print(vall)
            else:
                # calculate mae
                maelist = []
                for v1 in vall:
                    maes = []
                    for v2 in vall:
                        if metric == "mae":
                            mae = mean_absolute_error(v1, v2)
                        elif metric == "max":
                            v1 = np.array(v1)
                            v2 = np.array(v2)
                            mae = np.max(np.abs(v1-v2))
                        elif metric == "rmax":
                            v1 = np.array(v1)
                            v2 = np.array(v2)
                            mae = np.max(np.abs((v1-v2)/v2))
                        maes.append(mae)
                    maelist.append(maes)
                df_mae = pd.DataFrame(
                    maelist, index=dateprocesslist, columns=dateprocesslist)
                print(metric)
                print(df_mae)

            changed = True

    if changed:
        print(
            f"--->VALUE(s) CHANGED. go={go}, target={target}, subtarget={subtarget}")
    else:
        print(
            f"--->NO VALUE CHANGED. go={go}, target={target}, subtarget={subtarget}")
    print()
    return changed


@click.group()
@click.argument("action", type=click.Choice(['two', 'history']), default="history")
@click.option("--datelist", help="date list", show_default=True,
              default=["0830", "0902", "0905", "0909", "0915", "current"])
@click.option("--compiler", help="compiler", show_default=True, default="ifort")
@click.option("--prefix", help="Akaikkr path",  show_default=True, default="AkaiKKRprogram")
@click.option("--directory1", default=None)
@click.option("--directory2", default=None)
@click.option("--kkrtype", help="kkrtype",  default="akaikkr_cnd")
@click.option("--metric", help="metric",  show_default=True, default="max")
@click.option("--show_df", help="show_df",  show_default=True, type=bool, default=True)
@click.option("--resultjson", help="result file in the json format",
              show_default=True, default="result.json")
@click.pass_context
def cmd(ctx, action,
        datelist, compiler, prefix,
        directory1, directory2,
        kkrtype, metric, show_df, resultjson):
    """If action == history,
    you must specify directories by --datelist, --compiler, --prefix, --kkrtype.

    If action == two,
    you must specify directories by --directory1 --directory2, --kkrtype.
    """
    if action == "two":
        ctx.obj = {"directory1": directory1, "directory": directory2,
                   "kkrtype": kkrtype,
                   "metric": metric,
                   "show_df": show_df, "resultjson": resultjson}
        filepathlist, idlist = make_filepathlist_two(
            directory1, directory2, kkrtype, resultjson)
    elif action == "history":
        ctx.obj = {"datelist": datelist, "compiler": compiler, "prefix": prefix,
                   "kkrtype": kkrtype, "metric": metric,
                   "show_df": show_df, "resultjson": resultjson}
        filepathlist, idlist = make_filepathlist_history(
            prefix, datelist, compiler, kkrtype, resultjson)
    else:
        print("unknown action", action)
        raise ValueError("unknown action")
    ctx.obj.update({"filepathlist": filepathlist, "idlist": idlist})


def _go(ctx, go, target):
    if target is None:
        target_list = ["te", "tm", "spinlocalmoment", "orbitallocalmoment"]
    else:
        target_list = [target]
    subtarget = None
    notchanged_dic = OrderedDict()
    for target in target_list:
        changed = compare_them(go, target, subtarget, **ctx)
        notchanged_dic[target] = not changed
    notchanged_array = np.array(list(notchanged_dic.values()))
    if np.all(notchanged_array == True):
        print("ALL TEST PASSED")
    else:
        print("VALUE(s) CHANGED")
        print(notchanged_dic)


@cmd.command()
@click.option("--target", help="target property", default=None)
@click.pass_obj
def go(ctx, target):
    """read out_go.log"""
    go = "go"
    _go(ctx, go, target)


@cmd.command()
@click.option("--target", help="target property", default=None)
@click.pass_obj
def gofmg(ctx, target):
    """read gofmg.log"""
    go = "gofmg"
    _go(ctx, go, target)


@cmd.command()
@click.option("--target", help="target property", default=None)
@click.pass_obj
def fsm(ctx, target):
    """read gofmg.log"""
    go = "fsm"
    _go(ctx, go, target)


@cmd.command()
@click.option("--target", help="target property", default="Tc")
@click.pass_obj
def tc(ctx, target):
    """read out_tc.log"""
    go = "tc"
    subtarget = None
    compare_them(go, target, subtarget, **ctx)


@cmd.command()
@click.option("--target", help="target property", default=None)
@click.pass_obj
def jij(ctx, target):
    """read out_j30.log or out_j.log"""
    go = "jij"
    if target is None:
        target_list = ["Tc", "jij"]
    else:
        target_list = [target]
    subtarget = None
    for target in target_list:
        compare_them(go, target, subtarget, **ctx)


_DOS_CHOICE_ = ["up", "dn"]


@cmd.command()
@click.option("--subtarget", help="subtarget property",
              type=click.Choice(_DOS_CHOICE_),
              default=None)
@click.pass_obj
def dos(ctx, subtarget):
    """read out_dos.log"""
    go = "dos"
    notchanged_dic = OrderedDict()
    target_list = ["totalDOS", "pdos"]
    for target in target_list:
        if subtarget is None:
            subtarget_list = _DOS_CHOICE_
        else:
            subtarget_list = [subtarget]
        for subtarget in subtarget_list:
            changed = compare_them(go, target, subtarget, **ctx)
            notchanged_dic["_".join([target, subtarget])] = not changed

    notchanged_array = np.array(list(notchanged_dic.values()))
    if np.all(notchanged_array == True):
        print("ALL TEST PASSED")
    else:
        print("VALUE(s) CHANGED")
        print(notchanged_dic)


_AW_CHOICE_ = ["Aw_up", "Aw_dn"]


@ cmd.command()
@ click.option("--subtarget", help="subtarget property",
               type=click.Choice(_AW_CHOICE_),
               default=None)
@ click.pass_obj
def spc(ctx, subtarget):
    """read out_spc31.log"""
    go = "spc31"
    target = "Awk"
    if subtarget is None:
        subtarget_list = _AW_CHOICE_
    else:
        subtarget_list = [subtarget]
    notchanged_dic = OrderedDict()
    for subtarget in subtarget_list:
        changed = compare_them(go, target, subtarget, **ctx)
        notchanged_dic[subtarget] = not changed
    notchanged_array = np.array(list(notchanged_dic.values()))
    if np.all(notchanged_array == True):
        print("ALL TEST PASSED")
    else:
        print("VALUE(s) CHANGED")
        print(notchanged_dic)


if __name__ == "__main__":
    cmd()
