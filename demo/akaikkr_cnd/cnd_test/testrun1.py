# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import os
from akaikkr_cnd import AkaikkrJob
#
#
akaikkr_exe="../../specx"
#
# main function
def main():
	# fcc PdAg
	result=[]
	conc_list=[10,20,30,40,50,60,70,80,90]
	for conc in conc_list:
		directory="Pd"+str(100-conc)+"Ag"+str(conc)
		os.makedirs(directory,exist_ok=True)
		print("# fcc PdAd (non-magnetic, mjwasa)",directory)
		job=AkaikkrJob(directory)
		param=job.default
		param["brvtyp"]="fcc" ; param["a"]=1000000
		param["magtyp"]="nmag" ; param["sdftyp"]="mjwasa" ; param["record"]="init"
		param["edelt"]=1e-3 ; param["bzqlty"]=10
		param["ntyp"]=1 ; param["type"]=["PdAg"]
		param["ncmp"]=[2] ; param["rmt"]=[0.0]
		param["field"]=[0.0] ; param["mxl"]=[3]
		param["d"]=[[[0.0,0.0,0.0],[0.0,0.0,0.0]]]
		param["anclr"]=[[46,47]] ; param["conc"]=[[100-conc,conc]]
		param["natm"]=1
		param["atmicx"]=[
                    ["0.0a","0.0b","0.0c","PdAg"]
                    ]
		# go cal.
		job.make_inputcard(param,param["d"],"inputcard_go")
		job.run(akaikkr_exe,"inputcard_go","out_go.log")
		job.show_summary_go("out_go.log")
		# cnd cal.
		param["go"]=" cnd" ; param["record"]="2nd"
		param["edelt"]=1e-6 ; param["ewidth"]=0.01
		param["bzqlty"]=40
		job.make_inputcard(param,param["d"],"inputcard_cnd")
		job.run(akaikkr_exe,"inputcard_cnd","out_cnd.log")
		job.show_summary_cnd("out_cnd.log")
		result.append(job.get_result_testrun("out_cnd.log"))
		#
		#
	result_Pd90Ag10_ref={"res":10.1973}
	result_Pd80Ag20_ref={"res":17.41109}
	result_Pd70Ag30_ref={"res":24.02682}
	result_Pd60Ag40_ref={"res":26.69461}
	result_Pd50Ag50_ref={"res":21.47605}
	result_Pd40Ag60_ref={"res":15.5203}
	result_Pd30Ag70_ref={"res":12.92602}
	result_Pd20Ag80_ref={"res":9.12832}
	result_Pd10Ag90_ref={"res":4.82679}
	#
	#
	print("# resistivity difference between present and reference cals.")
	print("  Pd90Ag10:  {:.4e}".format(result[0]["res"]-result_Pd90Ag10_ref["res"]))
	print("  Pd80Ag20:  {:.4e}".format(result[1]["res"]-result_Pd80Ag20_ref["res"]))
	print("  Pd70Ag30:  {:.4e}".format(result[2]["res"]-result_Pd70Ag30_ref["res"]))
	print("  Pd60Ag40:  {:.4e}".format(result[3]["res"]-result_Pd60Ag40_ref["res"]))
	print("  Pd50Ag50:  {:.4e}".format(result[4]["res"]-result_Pd50Ag50_ref["res"]))
	print("  Pd40Ag60:  {:.4e}".format(result[5]["res"]-result_Pd40Ag60_ref["res"]))
	print("  Pd30Ag70:  {:.4e}".format(result[6]["res"]-result_Pd30Ag70_ref["res"]))
	print("  Pd20Ag80:  {:.4e}".format(result[7]["res"]-result_Pd20Ag80_ref["res"]))
	print("  Pd10Ag90:  {:.4e}".format(result[8]["res"]-result_Pd10Ag90_ref["res"]))
	#
	#
main()

