# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import os
import subprocess
import numpy as np
from pymatgen.core import Element
from akaikkr_cnd import AkaikkrJob
#
#
akaikkr_exe="../../../specx"
#
# main function
def main():
	# target list
	target_list=["Cu"]
	temp_list=[100,200,300,400,500,600,700,800,900,1000]
	#
	# open debye temp. list and make dic
	element_list=[]
	debye_temp_list=[]
	with open("debye_temp_list.txt") as f:
		next(f)
		data=[i.rstrip() for i in f.readlines()]
	for line in data:
		element_list.append(line.split()[1])
		debye_temp_list.append(float(line.split()[2])) # debye temp for 0K
		#debye_temp_list.append(float(line.split()[3])) # debye temp for 298K
	debye_temp_dic=dict(zip(element_list,debye_temp_list))
	#
	# open lattice const. list and make dic
	with open("lattice_constant_list.txt") as f:
		next(f)
		data=[i.rstrip() for i in f.readlines()]
	element_list=[]
	structure_list=[]
	lattice_const_list=[]
	for line in data:
		element_list.append(line.split()[0])
		structure_list.append(line.split()[1])
		lattice_const_list.append([float(line.split()[2]),float(line.split()[3])])
	structure_dic=dict(zip(element_list,structure_list))
	lattice_const_dic=dict(zip(element_list,lattice_const_list))
	#
	# target loop
	for target in target_list:
		elm=Element(target)
		z=elm.Z
		m=float(elm.atomic_mass)
		print("#",target)
		print(" ",structure_dic[target],lattice_const_dic[target],"(ang.)")
		#
    # debye temp.
		if (np.isnan(debye_temp_dic[target])):
			print("no debye temp. data",debye_temp_dic[target],target)
			sys.exit()
		debye_temp=debye_temp_dic[target]
		print("  atomic mass: ",m," debye temp: ",debye_temp)
		nrun=0
		#
		# temperature loop
		for temp in temp_list:
			nrun+=1
			directory=target+"/"+str(temp)
			os.makedirs(directory,exist_ok=True)
			job=AkaikkrJob(directory)
			param=job.default
			#
			param["a"]=lattice_const_dic[target][0]/0.529177
			param["c/a"]=lattice_const_dic[target][1]/lattice_const_dic[target][0]
			param["brvtyp"]=structure_dic[target]
			param["magtyp"]="nmag"
			param["sdftyp"]="mjwasa"
			param["record"]="init"
			param["edelt"]=1e-3
			param["ewidth"]=1.0
			param["bzqlty"]=10
			param["pmix"]=0.02
			param["ntyp"]=1
			param["type"]=[target]
			param["ncmp"]=[14]
			param["rmt"]=[0.0]
			param["field"]=[0.0]
			param["mxl"]=[4]
			param["anclr"]=[[z]*14] 
			param["conc"]=[[100]*14]
			if param["brvtyp"] == "hcp":
				param["natm"]=2
				param["atmicx"]=[
                        ["0.0a","0.0b","0.0c",target],
                        ["1/3a","2/3b","1/2c",target]
                        ]
			else:
				param["natm"]=1
				param["atmicx"]=[
                        ["0.0a","0.0b","0.0c",target]
                        ]
			#
			# mean square displacement by debye apporximation
			displace=job.avr_debye(m,debye_temp,temp)
			d=[job.make_displace(param,displace,14)]
			print("  temperature: ",temp,"(K)")
			print("  mean square displaement: ",displace,"(bohr)")
			print("  directory: ",directory)
			#
			# converged potential
			if (nrun != 1):
				cmd="cp {} {}".format(prev_dir+"/"+param["potentialfile"],directory+"/"+param["potentialfile"])
				subprocess.call(cmd,shell=True)

			# rough scf cal. for 1st concentrtion
			if (nrun == 1):
				job.make_inputcard(param,d,"inputcard_go")
				job.run(akaikkr_exe,"inputcard_go","out_go.log")

			# tight scf cal.
			param["record"]="2nd"
			param["bzqlty"]=20
			param["edelt"]=1e-4
			job.make_inputcard(param,d,"inputcard_go")
			job.run(akaikkr_exe,"inputcard_go","out_go.log")
      #
      # dos cal.
			param["go"]="dos"
			param["record"]="2nd"
			param["bzqlty"]=20
			param["edelt"]=1e-4
			job.make_inputcard(param,d,"inputcard_dos")
			job.run(akaikkr_exe,"inputcard_dos","out_dos.log")
			#
      # dos cal.
			param["go"]=" cnd"
			param["record"]="2nd"
			param["bzqlty"]=60
			param["edelt"]=1e-8
			param["ewidth"]=0.01
			job.make_inputcard(param,d,"inputcard_cnd")
			job.run(akaikkr_exe,"inputcard_cnd","out_cnd.log")

			prev_dir=directory
		

main()

