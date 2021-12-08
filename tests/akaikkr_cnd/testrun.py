# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import os
import sys
#cwd = os.getcwd()
#_path = os.path.join(os.getcwd(), "../../util")
#print("add python path",_path)
#sys.path.append(_path)
if "pyakaikkr" not in sys.modules:
    from pyakaikkr import *
if "akaikkr_testscript" not in sys.modules:
    from akaikkr_testscript import *


def make_exe():
    exe_dic = {}

    exe_list = []
    exe_list += [Cu_go, Cu_dos, Cu_spc ]
    exe_list += [Co_go, Co_fsm, Co_tc, Co_j30,  Co_dos, Co_spc  ]
    exe_list += [Fe_go, Fe_fsm, Fe_tc, Fe_j30,Fe_dos, Fe_spc ]
    exe_list += [Ni_go, Ni_fsm, Ni_tc, Ni_j30, Ni_dos, Ni_spc]
    exe_list += [NiFe_go, NiFe_fsm,NiFe_tc,NiFe_j30,  NiFe_dos, NiFe_cnd, NiFe_spc ]
    exe_list += [FeRh05Pt05_go, FeRh05Pt05_j30, FeRh05Pt05_fsm, FeRh05Pt05_tc, FeRh05Pt05_dos, FeRh05Pt05_spc, FeRh05Pt05_cnd]
    exe_list += [AlMnFeCo_bcc_go, AlMnFeCo_bcc_fsm, AlMnFeCo_bcc_gofmg, AlMnFeCo_bcc_tc, 
        AlMnFeCo_bcc_j30, 
        AlMnFeCo_bcc_dos, AlMnFeCo_bcc_cnd, AlMnFeCo_bcc_spc]
    exe_list += [Fe_lmd_go,  Fe_lmd_dos, Fe_lmd_spc]
    exe_list += [GaAs_go,  GaAs_dos, GaAs_spc ]
    exe_list += [Co2MnSi_go, Co2MnSi_fsm, Co2MnSi_tc, Co2MnSi_j30,
                  Co2MnSi_dos, Co2MnSi_spc ]
    exe_list += [SmCo5_oc_go, SmCo5_oc_fsm,  SmCo5_oc_tc, 
                SmCo5_oc_j30, SmCo5_oc_dos, SmCo5_oc_spc ]
    exe_dic.update({"all":exe_list})

    exe_list = []
    exe_list += [ Co_go, Co_j30 ]
    exe_dic.update({"Co": exe_list})

    exe_list = []
    exe_list += [ Fe_go, Fe_j30 ]
    exe_dic.update({"Fe": exe_list})

    exe_list = []
    exe_list += [ Ni_go, Ni_j30 ]
    exe_dic.update({"Ni": exe_list})

    exe_list = []
    exe_list += [ NiFe_go, NiFe_j30 ]
    exe_dic.update({"NiFe": exe_list})

    exe_list = []
    exe_list += [FeRh05Pt05_go, FeRh05Pt05_j30]
    exe_dic.update({"FeRhPt": exe_list})

    exe_list = []
    exe_list += [SmCo5_oc_go, SmCo5_oc_j30]
    exe_dic.update({"SmCo5_oc": exe_list})

    exe_list = []
    exe_list += [Cu_go, Cu_dos, Cu_spc ]
    exe_dic.update({"Cu": exe_list})

    exe_list = []
    exe_list += [Fe_lmd_go]
    exe_dic.update({"Fe_lmd": exe_list})

    exe_list = []
    exe_list += [Co2MnSi_go, Co2MnSi_j30]
    exe_dic.update({"Co2MnSi": exe_list})

    exe_list = []
    exe_list += [Cu_go, Cu_dos, Cu_spc ]
    exe_list += [Co_go, Co_fsm, Co_tc, Co_j30,  Co_dos, Co_spc  ]
    exe_list += [Fe_go, Fe_fsm, Fe_tc, Fe_j30,Fe_dos, Fe_spc ]
    exe_list += [Ni_go, Ni_fsm, Ni_tc, Ni_j30, Ni_dos, Ni_spc]
    exe_dic.update({"2": exe_list})

    return exe_dic

if __name__ == "__main__":
    akaikkr_exe = "specx"
    fmg_exe = "fmg"
    exe_dic = make_exe()
    all_go(akaikkr_exe, fmg_exe=fmg_exe, exe_dic=exe_dic, displc=True)

