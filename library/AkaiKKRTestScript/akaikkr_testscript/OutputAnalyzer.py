# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


from copy import deepcopy
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error

_SUCCESS_ = "O"
_FAILED_ = "X"
_NO_REF_ = "-"

_FILL_STR_ = "-"
_PAD_STR_ = "  "

_CURRENT_ = "current"
_REFERENCE_ = "reference"

_UPDN_LIST_ = ["up", "dn"]

_UNNUMBERED_TARGET_ = ["te", "tm", "Tc", "resis"]
_NUMBERED_TARGET_ = ["spinlocalmoment", "orbitallocalmoment",
                     "cnd"]


def show_chk_legend():
    print("{}: passed. {}: failed, {}: no reference".format(
        _SUCCESS_, _FAILED_, _NO_REF_))
    print()


def chk_dic_success_all(chk_dic):
    chk_array = np.array([v == _SUCCESS_ for v in chk_dic.values()])
    return chk_array


def _go_fix_localmoment(result, ref, shortname=False):
    """ref can be None"""
    result = deepcopy(result)
    ref = deepcopy(ref)
    df1 = pd.DataFrame([[x for x in result.values()]], index=[_CURRENT_],
                       columns=result.keys())
    if ref is not None:
        df2 = pd.DataFrame([[x for x in result.values()]], index=[_REFERENCE_],
                           columns=result.keys())
    else:
        df2 = None
    return df1, df2
    if True:
        for obj in [result, ref]:
            if obj is None:
                continue
            for label in _NUMBERED_TARGET_:
                if label in obj:
                    v = obj[label]
                    if isinstance(v, list):
                        del obj[label]
                        for i, r in enumerate(v):
                            name = "{}{}".format(label, i+1)
                            obj[name] = r

    result.update({"target": _CURRENT_})
    if ref is not None:
        ref.update({"target": _REFERENCE_})
    df1 = pd.DataFrame(result, dtype=object, index=[0])
    try:
        df2 = pd.DataFrame(ref, dtype=object)
    except ValueError:
        df2 = pd.DataFrame(ref, dtype=object, index=[0])

    # replace long name with short name
    if shortname:
        for df in [df1, df2]:
            _make_df_shortname(df)
    return df1, df2


class ThresBase:
    def __init__(self, thres):
        self.thres = thres
        thres_key_type = []
        thres_key_name = []
        for thres_key, thres_value in thres.items():
            s = thres_key.split("_")
            thres_key_type.append(deepcopy(s[0]))
            thres_key_name.append(deepcopy(s[1]))
        thres_key_type = list(set(thres_key_type))
        thres_key_name = list(set(thres_key_name))
        thres_key_type.sort()
        thres_key_name.sort()
        self.thres_key_type = thres_key_type
        self.thres_key_name = thres_key_name

    @property
    def key_type(self):
        return self.thres_key_type

    @property
    def key_name(self):
        return self.thres_key_name


def _make_df_diff(key, result, ref, thres):
    df1, df2 = _go_fix_localmoment(result, ref)
    if ref is None:
        return df1
    df = pd.concat([df1, df2])
    res = {}
    thresbase = ThresBase(thres)
    for thres_key_type in thresbase.key_type:
        res[thres_key_type] = {}
        thres_key_type_th = thres_key_type+"_th"
        res[thres_key_type_th] = {}

    for thres_key, thres_value in thres.items():
        s = thres_key.split("_")
        thres_key_type = deepcopy(s[0])
        thres_key_name = deepcopy(s[1])
        thres_key_type_th = thres_key_type+"_th"

        for name in df.columns:
            if thres_key_name == name:
                v_current = df.loc[_CURRENT_, name]
                if isinstance(v_current, list):
                    v_current = np.array(v_current)
                v_reference = df.loc[_REFERENCE_, name]
                if isinstance(v_reference, list):
                    v_reference = np.array(v_reference)
                if thres_key_type == "diff":
                    value = np.abs(v_current - v_reference)
                    res[thres_key_type].update({name: value})
                    res[thres_key_type_th].update({name: thres_value})
                elif thres_key_type == "rdiff":
                    value = np.abs(v_current - v_reference) /\
                        np.abs(v_reference)
                    res[thres_key_type].update({name: value})
                    res[thres_key_type_th].update({name: thres_value})

                elif thres_key_type == "and":
                    value = v_current and v_reference
                    res[thres_key_type].update({name: value})
                    res[thres_key_type_th].update({name: thres_value})

                else:
                    print("unkown type", type)
                    print("thres", thres)
                    raise ValueError
    df_list = []
    for res_key, res_value in res.items():
        _df = pd.DataFrame(
            [[v for v in res_value.values()]], columns=res_value.keys(), index=[res_key])
        df_list.append(_df)
    df_comp = pd.concat(df_list, axis=0)
    df_all = pd.concat([df, df_comp], axis=0)
    return df_all


def _make_df_shortname(df):
    return df
    df = df.copy()
    col2 = []
    col = list(df.columns)
    for s in col:
        for lmtarget, abbrev in zip(["spinlocalmoment", "orbitallocalmoment", "threads"],
                                    ["sl", "ol", "thrd"]):
            if lmtarget in s:
                s = s.replace(lmtarget, abbrev)
        col2.append(s)
    df.columns = col2
    return df


def __df_str_add_spc(df, pad=_PAD_STR_):
    lines = []
    for x in df.__str__().splitlines():
        lines.append(pad+x)
    return "\n".join(lines)


class DiffResult:
    def __init__(self, key, result, ref, thres):
        self.key = key
        self.result = result
        self.ref = ref
        self.thres = thres
        df = _make_df_diff(key, result, ref, thres)
        self.df = df

    def process(self, ):
        df = self.df
        thres = self.thres
        if _REFERENCE_ not in df.index:
            return df, None

        res_dic = {}
        for thres_op_label in thres.keys():
            s = thres_op_label.split("_")
            op = s[0]
            label = s[1]
            op_th = op+"_th"
            value = df.loc[op, label]
            if isinstance(value, list):
                value = np.array(value)
            thvalue = df.loc[op_th, label]
            if op == "and":
                flag = value == thvalue
            else:
                flag = value < thvalue
            if isinstance(flag, np.ndarray):
                flag = np.all(flag == True)
            res_dic[thres_op_label] = flag
        df_chk = pd.DataFrame(res_dic, index=["chk"])

        return df, df_chk


class DiffVector:
    def __init__(self, df, thres):
        self.df = df
        self.thres = thres
        self. difflabel = "diff"

    def process(self):
        """make dataframe
        """
        df = self.df
        difflabel = self.difflabel
        v1 = df.loc[_CURRENT_, :].astype(float).values
        v2 = df.loc[_REFERENCE_, :].astype(float).values
        diffvalue = np.abs(v1-v2).tolist()
        _df = pd.DataFrame([v1.tolist(), v2.tolist(), diffvalue], columns=df.columns,
                           index=[_CURRENT_, _REFERENCE_, difflabel])

        self.df = _df  # update self.df
        return _df

    def evaluate_distance(self):
        _df = self.df
        difflabel = self.difflabel
        thres = self.thres

        diff_dic = {}
        for type_ in thres.keys():
            if type_ == "diff_max":
                value = np.max(_df.loc[difflabel, :].values)
            elif type_ == "mae":
                value = mean_absolute_error(_df.loc[_CURRENT_, :].values,
                                            _df.loc[_REFERENCE_, :].values)
            else:
                print("unknown thres.keys()", thres.keys())
                raise ValueError
            diff_dic.update({type_: value})
        df_diff = pd.DataFrame(
            diff_dic, dtype='object', index=["value"])

        th = {}
        for label in ["diff_max", "mae"]:
            th[label] = thres[label]
        df_th = pd.DataFrame(th, index=["thres"])
        df_diff = pd.concat([df_diff, df_th], axis=0)

        chk = {}
        for diffkey in ["diff_max", "mae"]:
            chk[diffkey] = df_diff.loc["value",
                                       diffkey] < df_diff.loc["thres", diffkey]
        df_chk = pd.DataFrame(chk, index=["chk"])

        df_diff = pd.concat([df_diff, df_chk], axis=0)
        return df_diff


def _sort_types_inside_row(df):
    """exchange type1 and type2 as type1<type2 to compare with another data.
    and add pair column.
    """
    df = df.copy()
    type1 = df["type1"].values
    type2 = df["type2"].values
    t1t2_list = []
    for t1, t2 in zip(type1, type2):
        _t1t2 = [t1, t2]
        _t1t2.sort()
        t1t2_list.append(_t1t2)
    df[["type1", "type2"]] = t1t2_list

    if "pair" in list(df.columns):
        del df["pair"]

    # add pair column
    comp1 = df["comp1"].values
    comp2 = df["comp2"].values
    type1 = df["type1"].values
    type2 = df["type2"].values
    typepair = []
    for t1, t2, c1, c2 in zip(type1, type2, comp1, comp2):
        typepair.append("-".join([t1, t2, c1, c2]))
    df_pair = pd.DataFrame({"pair": typepair})

    jijdf = pd.concat([df, df_pair], axis=1)
    return jijdf


def _make_jij_dataframe(result_jij, ref_jij, target="J_ij(meV)"):
    # make jij dataframe
    df_result_jij = pd.DataFrame(result_jij[1:], columns=result_jij[0])
    df_result_jij = df_result_jij[[
        "comp1", "comp2", "J_ij", "J_ij(meV)", "pair"]]

    if ref_jij is not None:
        df_ref_jij = pd.DataFrame(ref_jij[1:], columns=ref_jij[0])
        df_ref_jij = _sort_types_inside_row(df_ref_jij)

        # sort values to compare the result with another dataframe.
        df_ref_jij.sort_values(by="distance", inplace=True)
        df_ref_jij = df_ref_jij[[
            "comp1", "comp2", "J_ij", "J_ij(meV)", "pair"]]
        df_ref_jij.reset_index(drop=True)

        zipped_verson_df_list = zip([_CURRENT_, _REFERENCE_], [
            df_result_jij, df_ref_jij])
    else:
        zipped_verson_df_list = zip([_CURRENT_], [df_result_jij])

    # replace J_ij* with {_CURRENT_}_J_ij*.
    # replace distace with {_CURRENT_}_distance.
    for version, _df in zipped_verson_df_list:
        col2 = {}
        for col in _df.columns:
            if col.startswith("J_ij"):
                col2[col] = "{}_{}".format(version, col)
            if col.startswith("distance"):
                col2[col] = "{}_{}".format(version, col)
        _df.rename(columns=col2, inplace=True)

    if True:
        print()
        print("debug df_result_jij")
        print(df_result_jij)

        print("debug df_ref_jij")
        print(df_ref_jij)
        print()

    if ref_jij is not None:
        # merge
        col = ["comp1", "comp2", 'pair']
        df = df_result_jij.merge(df_ref_jij, on=col)
    else:
        df = df_result_jij

    # add pair-comp1-comp2 field
    paircomp = []
    for comp1, comp2, pair in zip(df["comp1"], df["comp2"], df["pair"]):
        paircomp.append("{}_{}_{}".format(pair, comp1, comp2))
    df["typecomp"] = paircomp

    if "{}_{}".format(_REFERENCE_, target) in df.columns:
        col = ["{}_{}".format(_CURRENT_, target), "{}_{}".format(
            _REFERENCE_, target), "typecomp"]
    else:
        col = ["{}_{}".format(_CURRENT_, target), "typecomp"]
    _df = df[col].set_index("typecomp").T

    # rename index
    col2 = {}
    for x in col:
        s = x.split("_")
        col2[x] = s[0]
    _df = _df.rename(index=col2)
    return _df


def _spactra_df(result_totaldos, ref_totaldos,
                updn_list=_UPDN_LIST_, updn_label_list=_UPDN_LIST_):

    if ref_totaldos is not None:
        zipped_version_dos = zip([_CURRENT_, _REFERENCE_],
                                 [result_totaldos, ref_totaldos])
    else:
        zipped_version_dos = zip([_CURRENT_, ],
                                 [result_totaldos, ])
    for version, obj in zipped_version_dos:
        for updn, updnlabel in zip(updn_list, updn_label_list):
            if updnlabel in obj:
                dos = obj[updnlabel]
                if len(dos) > 0:
                    name = "{}_{}".format(version, updn)
                    obj[name] = obj.pop(updnlabel)
                else:
                    del obj[updn]

    df_result_totaldos = pd.DataFrame(result_totaldos)
    if ref_totaldos is not None:
        df_ref_totaldos = pd.DataFrame(ref_totaldos)

        if "energy" in df_ref_totaldos:
            # confirm that their eneriges are the same
            energy_result = df_result_totaldos.loc[:, "energy"].values
            energy_ref = df_ref_totaldos.loc[:, "energy"].values
            if np.all(energy_result == energy_ref):
                df = df_result_totaldos.merge(df_ref_totaldos, on="energy")
            else:
                print(
                    "energies are different between the current calculation and reference.")
                print("no check applied.")
                df = df_result_totaldos
        else:
            print("no energy in reference data. but continue.")
            df = pd.concat([df_result_totaldos, df_ref_totaldos], axis=1)
    else:
        df = df_result_totaldos

    # df is now
    #      energy  _CURRENT__up  _CURRENT__dn  _REFERENCE__up  _REFERENCE__dn
    # 0    -1.495     0.00094     0.00089       0.00094       0.00089
    # 1    -1.485     0.00096     0.00090       0.00096       0.00090
    #
    # split them into up and dn
    df_dos_updn = {}
    for updn in updn_list:
        # search **_up" in columns
        col = []
        for s in df.columns:
            if s.endswith("_"+updn):
                col.append(s)

        if len(col) > 0:
            col.append("energy")
            _df = df.loc[:, col].set_index("energy").T
            index2 = {}
            for s in _df.index:
                s2 = s.replace("_"+updn, "")
                index2[s] = s2
            _df = _df.rename(index=index2)
            df_dos_updn[updn] = _df
    return df_dos_updn


def go_diff_msg(key, result, ref,
                thres={"and_conv": True,
                       "rdiff_te": 1e-11,
                       "diff_tm": 1e-5,
                       "diff_spinlocalmoment": 1e-5,
                       "diff_orbitallocalmoment": 1e-5, }):
    diffresult = DiffResult(key, result, ref, thres)
    df, df_chk = diffresult.process()
    # df.dropna(axis=0, inplace=True)
    df.fillna(_FILL_STR_, inplace=True)

    thresbase = ThresBase(thres)
    col = thresbase.key_name
    print(__df_str_add_spc(df[col]))
    print()
    if _REFERENCE_ not in df.index:
        return _NO_REF_

    print(_PAD_STR_, "summary")
    print(__df_str_add_spc(df_chk.T))
    print()

    chk = df_chk.loc["chk", :].values
    finalchk = np.all(chk == True)
    print(_PAD_STR_, "final:", finalchk)
    if ref is not None:
        if finalchk:
            return _SUCCESS_
        else:
            return _FAILED_
    else:
        return _NO_REF_


def cnd_diff_msg(key, result, ref, thres={"rdiff_resis": 1e-5,
                                          "rdiff_cnd": 1e-5}):
    return go_diff_msg(key, result, ref, thres=thres)


def tc_diff_msg(key, result, ref, thres={"rdiff_Tc": 1e-5}):
    return go_diff_msg(key, result, ref, thres=thres)


def j_diff_msg(key, result, ref,
               thres={"rdiff_Tc": 1e-5, },
               thres_vector={"diff_max": 1e-4, "mae": 1e-4, }):

    target = "jij"
    vector_df_title = "J_ij(meV) at the shortest distance"
    result_jij = result.pop(target, None)
    if ref is not None:
        ref_jij = ref.pop(target, None)
    else:
        ref_jij = None
    # please refer ref_jij, not ref from now on.

    if thres is not None:
        chk1_msg = go_diff_msg(key, result, ref, thres=thres)
        if chk1_msg == _SUCCESS_:
            chk1 = True
        elif chk1_msg == _FAILED_:
            chk1 = False
        elif chk1_msg == _NO_REF_:
            chk1 = True
        else:
            print("unknown return value, chk1", chk1_msg)
            raise ValueError
    else:
        chk1 = True

    df_jij = _make_jij_dataframe(result_jij, ref_jij)

    if _REFERENCE_ in df_jij.index:
        diffvector = DiffVector(df_jij, thres_vector)
        df_jij = diffvector.process()
        print(_PAD_STR_, vector_df_title)
        print(__df_str_add_spc(df_jij))
        print()
    else:
        print(_PAD_STR_, vector_df_title)
        print(__df_str_add_spc(df_jij))
        print()
        return _NO_REF_

    df_jij_chk = diffvector.evaluate_distance()
    print(_PAD_STR_, "{} summary".format(target))
    print(__df_str_add_spc(df_jij_chk))
    print()

    # must add df_tc_chk
    chk2 = np.all(df_jij_chk.loc["chk", :].values == True)
    chk = np.array([chk1, chk2])
    finalchk = np.all(chk == True)
    print(_PAD_STR_, "final:", finalchk)
    if finalchk:
        return _SUCCESS_
    else:
        return _FAILED_


def vector_diff_msg(result_totaldos, ref_totaldos,
                    vector_df_title, thres_vector, updn_label_list=_UPDN_LIST_,):

    df_dos_updn = _spactra_df(
        result_totaldos, ref_totaldos, updn_label_list=updn_label_list)
    if _REFERENCE_ in df_dos_updn["up"].index:
        df_chk_updn = {}
        for (updn, df_dos), vector_df_title_spin in zip(df_dos_updn.items(), vector_df_title):
            diffvector = DiffVector(df_dos, thres_vector)
            df_dos = diffvector.process()
            print(_PAD_STR_, vector_df_title_spin)
            print(__df_str_add_spc(df_dos))
            print()
            if _REFERENCE_ in df_dos.index:
                df_chk = diffvector.evaluate_distance()
                print(__df_str_add_spc(df_chk))
                print()
                df_chk_updn[updn] = df_chk
    else:
        for (updn, df_dos), vector_df_title_spin in zip(df_dos_updn.items(), vector_df_title):
            print(_PAD_STR_, vector_df_title_spin)
            print(__df_str_add_spc(df_dos))
            print()
        return _NO_REF_

    chk_list = []
    for df_chk in df_chk_updn.values():
        chk = df_chk.loc["chk", :].values
        chk_list.append(np.all(chk == True))
    chk_list = np.array(chk_list)
    finalchk = np.all(chk_list == True)
    print(_PAD_STR_, "vector diff:", finalchk)
    print()
    if finalchk:
        return _SUCCESS_
    else:
        return _FAILED_


def dos_diff_msg(key, result, ref, target_list=["totalDOS", "pdos"],
                 updn_label_list=_UPDN_LIST_,
                 thres_vector={"diff_max": 2.5e-3,
                               "mae": 2.5e-4, }):
    result = deepcopy(result)

    chk_list = []
    for target in target_list:
        result_totaldos = result.pop(target)
        if ref is not None:
            ref = deepcopy(ref)
            ref_totaldos = ref.pop(target, None)
        else:
            ref_totaldos = None

        if target == "totalDOS":
            vector_df_title = []
            for updn in updn_label_list:
                vector_df_title.append("{} ({})".format(target, updn))
        elif target == "pdos":
            vector_df_title = []
            atom = result_totaldos["atom"]
            l = result_totaldos["l"]
            for updn in updn_label_list:
                vector_df_title.append(
                    "{} ({}), atom={}, l={}".format(target, updn, atom, l))

        chk = vector_diff_msg(result_totaldos, ref_totaldos,
                              vector_df_title=vector_df_title, thres_vector=thres_vector)
        if chk == _SUCCESS_ or chk == _NO_REF_:
            chk = True
        else:
            chk = False
        chk_list.append(chk)

    chk_list = np.array(chk_list)
    finalchk = np.all(chk_list == True)
    print(_PAD_STR_, "final:", finalchk)
    print()
    if finalchk:
        return _SUCCESS_
    else:
        return _FAILED_


def aw_diff_msg(key, result, ref,
                target="Awk",
                updn_label_list=["Aw_up", "Aw_dn"],
                thres_vector={"diff_max": 3e-3,
                              "mae": 1e-4, }):

    result = deepcopy(result)
    resultAwk = result.pop(target, None)
    if ref is not None:
        ref = deepcopy(ref)
        refAwk = ref.pop(target, None)
    else:
        refAwk = None

    updn_label_list = []
    for updn in _UPDN_LIST_:
        updn_label_list.append("Aw_{}".format(updn))

    vector_df_title = []
    for updn in _UPDN_LIST_:
        vector_df_title.append("Aw(w,k) k={}, {}".format(
            resultAwk["Awk_kpoint"], updn))

    return vector_diff_msg(resultAwk, refAwk, updn_label_list=updn_label_list,
                           vector_df_title=vector_df_title, thres_vector=thres_vector)
