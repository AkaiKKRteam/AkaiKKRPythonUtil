# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

import sys
import subprocess
import numpy as np
from scipy import integrate

class AkaikkrJob:
	#
	#
	def __init__(self,path_dir):
		# path of running directory
		self.path_dir=path_dir
		# default parameters for AkaiKKR
		self.default={
    "go": "go", "potentialfile": "pot.dat",
    "brvtyp": "bcc", "a": 0, "c/a": 1.0, "b/a": 1.0,
    "alpha": 90, "beta": 90, "gamma": 90,
    "edelt": 1e-3, "ewidth": 1.0,
    "reltyp": "sra", "sdftyp": "mjw", "magtyp": "mag", "record": "2nd",
    "outtyp": "update", "bzqlty": 6, "maxitr": 200, "pmix": "0.02",
    "ntyp": 1,
    "rmt": [1.0],
    "field": [0.0],
    "mxl": [3],
    "type": ["Fe"],
    "ncmp": [1],
    "d": [[[0.0,0.0,0.0]]],
    "anclr": [[26]],
    "conc": [[100]],
    "natm": 1,
    "atmicx": [
              ["0.00a","0.00b","0.00c","Fe"],
              ],
    }
	#
	#
	def make_inputcard(self,dic,d,inputcard):
		# make inputcard from param
		#
		def from_keylist(dic,keylist):
			s=[] ; result=[]
			head=["#---"]
			for key in keylist:
				head.append(key)
				s.append(str(dic[key]))
			result.append(" ".join(head))
			result.append(" ".join(s))
			return result
		#
		def make_rbasis(dic):
			result=[]
			for vec in ["r1","r2","r3"]:
				s=dic[vec]
				result.append(" ".join([str(i) for i in s]))
			result.append(str(dic["a"]))
			return result
		#
		def make_type_card(dic):
			result=[]; nn=0
			for i,j,k,l,m, in zip(dic["rmt"],dic["field"],dic["mxl"],dic["type"],dic["ncmp"]):
				result.append("#--- type ncmp")
				result.append(" ".join([" ",l,str(m)]))
				result.append("#- rmt field mxl")
				result.append(" ".join([str(i),str(j),str(k)]))
				result.append("#- anclr conc")
				for n in range(dic["ncmp"][nn]):
					result.append(" ".join([" ".join(map(str,d[nn][n])),str(dic["anclr"][nn][n]),str(dic["conc"][nn][n])]))
					#result.append(" ".join(["0.0 0.0 0.0",str(dic["anclr"][n]),str(dic["conc"][n])]))
					#result.append(" ".join([str(dic["anclr"][nn][n]),str(dic["conc"][nn][n])]))
				nn+=1
			return result
		#
		def make_atom_card(dic):
			result=[]
			result.append("#--- atmicx atmtyp")
			for n in range(dic["natm"]):
				#s=[" ".join(dic["atmicx"][n]),dic["atmtyp"][n]]
				s=dic["atmicx"][n]
				result.append(" ".join(s))
			return result
		#
		card=[]
		card+=from_keylist(dic,["go","potentialfile"])
		if dic["brvtyp"]=="aux":
			card+=from_keylist(dic,["brvtyp"])
			card+=make_rbasis(dic)
		else:
			card+=from_keylist(dic,["brvtyp","a","c/a","b/a","alpha","beta","gamma"])
		card+=from_keylist(dic,["edelt","ewidth","reltyp","sdftyp","magtyp","record"])
		card+=from_keylist(dic,["outtyp","bzqlty","maxitr","pmix"])
		card+=from_keylist(dic,["ntyp"])
		card+=make_type_card(dic)
		card+=from_keylist(dic,["natm"])
		card+=make_atom_card(dic)
		with open(self.path_dir+"/"+inputcard, mode="w") as f:
			f.write("\n".join(card))
	#
	#
	def run(self,akaikkr_exe,infile,outfile):
		# execute akaikkr
		cmd = "cd {}; {} < {} > {}".format(self.path_dir,akaikkr_exe,infile,outfile)
		subprocess.call(cmd,shell=True)
	#
	#
	def copy_potential(self,potentialfile1,potentialfile2):
		# copy potential file: file1 => file2
		cmd = "cd {}; cp {} {}".format(self.path_dir,potentialfile1,potentialfile2)
		subprocess.call(cmd,shell=True)
	#
	#
	def delete_potential(self,potentialfile):
		# delete potential file
		cmd = "cd {}; rm -rf {}".format(self.path_dir,potentialfile)
		subprocess.call(cmd,shell=True)
	#
	#
	def show_summary_go(self,outfile):
		# get summary of go calculation
		print("  directory: ",self.path_dir,"  output file: ",outfile)
		print("  converged: ",self.check_convergence_go(outfile))
		print("  total energy: ",self.get_total_energy(outfile), " Ry")
		print("  total moment: ",self.get_total_moment(outfile), " muB")
		print("  local moment: ",self.get_local_moment(outfile), " muB\n")
	#
	#
	def show_summary_cnd(self,outfile):
		# get summary of go calculation
		print("  directory: ",self.path_dir,"  output file: ",outfile)
		print("  resistivity: ",self.get_resistivity(outfile), " (micro ohm cm)\n")
	#
	#
	def get_result_testrun(self,outfile):
		# get result for testrun.
		# return dictionary
		dic={"te": self.get_total_energy(outfile),"rms": self.get_rms_error(outfile)[-1],
         "tm": self.get_total_moment(outfile),"threads": self.get_threads_openmp(outfile),
         "res": self.get_resistivity(outfile)}
		return dic
	#
	#
	def check_convergence_go(self,outfile):
		# check convergence for go calcualtion
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		converged=False
		for line in data:
			if "*** no convergence" in line: break
			if "sbtime report" in line: converged=True
		return converged
	#
	#
	def get_threads_openmp(self,outfile):
		# get the number of threads for OpenMP
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "threads" in line:
				threads=int(line.split("(")[1].split()[0])
		return threads
	#
	#
	def get_rms_error(self,outfile):
		# get information of rms-error
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		rms_error=[]
		for line in data:
			if "te=" in line:
				rms_error.append(float(line.split("err=")[1].split()[0]))
		return rms_error
	#
	#
	def get_lattice_constant(self,outfile):
		# get lattice constant
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "bravais=" in line:
				lattice_constant=float(line.split("a=")[1].split()[0])
		return lattice_constant
	#
	#
	def get_unitcell_volume(self,outfile):
		# get volume
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "unit cell volume=" in line:
				unitcell_volume=float(line.rstrip("\n").split("volume=")[1][:-6])
		return unitcell_volume
	#
	#
	def get_ewidth(self,outfile):
		# get ewidth
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "ewidth=" in line:
				ewidth=float(line.split("ewidth=")[1].split()[0])
		return ewidth
	#
	#
	def get_edelt(self,outfile):
		# get edelt
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "edelt=" in line:
				edelt=float(line.split("edelt=")[1].split()[0])
		return edelt
	#
	#
	def get_fermi_level(self,outfile):
		# get Fermi level
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "ef=" in line:
				ef_up=float(line.split()[1])
				ef_dn=float(line.split()[2])
				ef=(ef_up+ef_dn)/2
				break
		return ef
	#
	#
	def get_total_energy(self,outfile):
		# get total energy
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "total energy=" in line:
				te=float(line.split("total energy=")[-1].rstrip())
				break
		return te
	#	
	#
	def get_total_moment(self,outfile):
		# get total moment
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		h_moment=[]
		for line in data:
			if "te=" in line:
				h_moment.append(float(line.split("moment=")[1].split()[0]))
		return h_moment[-1]
	#
	#
	def get_local_moment(self,outfile,mode="spin"):
		# get information of local moment(spin or orbital)
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		spin_moment=[]
		orbital_moment=[]
		for line in data:
			if "spin moment=" in line:
				spin_moment.append(float(line.split("spin moment=")[1].split()[0]))
				orbital_moment.append(float(line.split("orbital moment=")[1].split()[0]))
		if mode == "spin":
			return spin_moment
		else:
			return orbital_moment
	#
	#
	def get_type_charge(self,outfile):
		# get total charge of type
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		charge=[]
		for line in data:
			if "total charge=" in line:
				charge.append(float(line.split("total charge=")[1].split()[0]))
		return charge
	#
	#
	def get_curie_temperature(self,outfile):
		# get curie temperature
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "Tc (in mean field approximation) =" in line:
				tc=float(line.split("Tc (in mean field approximation) =")[-1].rstrip()[:-1].strip())
				break
		return tc
	#
	#
	def get_resistivity(self,outfile):
		# get resistivity
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for line in data:
			if "resistivity" in line:
				resistivity=float(line.split()[1])
				break
		return resistivity 
	#
	#
	def get_conductivity_spin(self,outfile):
		# get resistivity
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		conductivity_spin=[]
		for line in data:
			if "cnd1" in line:
				conductivity_spin.append(float(line.split()[1]))
			if "cnd2" in line:
				conductivity_spin.append(float(line.split()[1]))
				break
		return conductivity_spin
	#
	#
	def check_core_level(self,outfile):
		# check core states
		core_state=["3d","4d","4f"]
		core_exist=[] ; core_level=[]
		with open(self.path_dir+"/"+outfile) as f:
			data=f.readlines()
		for core in core_state:
			dummy1=False ; dummy2=[]
			for line in data:
				if "Ry("+core+")" in line:
					dummy1=True
					dummy2.append(float(line.split("Ry("+core+")")[0].split()[-1]))
			core_exist.append(dummy1)
			core_level.append(dummy2)
		return core_exist,core_level
	#
	#
	def cut_dos_beta(self,dosfile,keyword,save=False):
		# cut dos (beta version)
		with open(self.path_dir+"/"+dosfile) as f:
			data=f.readlines()
		e=[]
		dos_up=[]
		dos_dn=[]
		nspin=0
		check=False
		for line in data:
			line=line.strip()
			if keyword in line:
				check=True ; nspin+=1
			if not line: check=False
			if check and keyword not in line:
				if nspin==1:
					e.append(float(line.split()[0]))
					dos_up.append(float(line.split()[1]))
				else:
					dos_dn.append(-float(line.split()[1]))
		if save:
			with open(self.path_dir+"/"+keyword.replace(" ","_")+".dat","a") as f:
				for i,j,k in zip(e,dos_up,dos_dn):
					f.write(" ".join([str(m) for m in [i,j,k]])+"\n")
		else:
			return e,dos_up,dos_dn
	#
	#
	def cut_dos(self,dosfile,keyword,save=False):
		# cut dos
		with open(self.path_dir+"/"+dosfile) as f:
			data=f.readlines()
		dos_up=[]
		dos_dn=[]
		nspin=0
		check=False
		for line in data:
			line=line.strip()
			if keyword in line:
				check=True ; nspin+=1
			if not line: check=False
			if check and keyword not in line:
				if nspin==1:
					dos_up.append(line)
				else:
					dos_dn.append(line)
		if save:
			with open(self.path_dir+"/"+keyword.replace(" ","_")+"_up.dat","w") as f:
				f.write("\n".join(dos_up)+"\n")
			with open(self.path_dir+"/"+keyword.replace(" ","_")+"_dn.dat","w") as f:
				f.write("\n".join(dos_dn)+"\n")
		else:
			return dos_up,dos_dn
	#
	#
	def cut_jij(self,jijfile,keyword,comp1=1,comp2=1,save=False):
		# cut jij
		with open(self.path_dir+"/"+jijfile) as f:
			data=f.readlines()
		jij=[]
		check=False
		for line in data:
			line=line.strip()
			if keyword in line: check=True
			if not line: check=False
			if check and keyword not in line and "index" not in line:
				if int(line.split()[3]) == comp1 and int(line.split()[4]) == comp2:
					jij.append(line)
		if save:
			outfile=self.path_dir+"/"+keyword+"_"+str(comp1)+"-"+str(comp2)+".dat"
			with open(outfile,"w") as f:
				f.write("\n".join(jij)+"\n")
		else:
			return jij
	#
	#
	def avr_debye(self,m,debye_temp,temp):
		# mean square dispalcement by debye approximation
		# define parameters
		kB=6.3336e-06            # Boltzmann constant, Ry/K
		me=9.10938356            # electron mass, kg (10**(-31))
		u_to_kg=1.660540         # conversion u to kg (10**(-27))
		#
		# convert atomic mass unit to atomic unit (i.e., me=0.5)
		m=m*u_to_kg/(2.0*me)*1.0e4
		#
		# define parameters for numerical integration
		mesh=10000               # mesh for integration
		xinit=1e-6               # initial point for integration
		#
		# debey function
		tau=debye_temp/temp
		t=np.linspace(xinit,tau,mesh)
		f=t/(np.exp(t)-1)
		D1=integrate.trapz(f,t)/tau
		#D1=integrate.simps(f,t)/tau
		#
		# mean square displacement
		fac=9.0/(m*kB*debye_temp)
		u2=fac*(D1/tau) # zero point energy is ignored
		#u2=fac*(D1/tau+0.25)
		u=np.sqrt(u2)
		return u
	#
	#
	def make_displace(self,dic,disp,ndir=14):
		u=disp/dic["a"]
		d2=[[0.0,0.0,u],[0.0,0.0,-u]]
		#d2=[[u,0.0,0.0],[-u,0.0,0.0]]
		d6=[[u,0.0,0.0],[-u,0.0,0.0],[0.0,u,0.0],[0.0,-u,0.0],[0.0,0.0,u],[0.0,0.0,-u]]
		u=disp/dic["a"]/np.sqrt(3)
		d8=[[u,u,u],[-u,-u,u],[u,-u,u],[-u,u,u],[u,u,-u],[-u,-u,-u],[u,-u,-u],[-u,u,-u]]
		if ndir == 2:
			return d2
		elif ndir == 6:
			return d6
		elif ndir == 8:
			return d8
		elif ndir == 14:
			return d6+d8
		else:
			print("bad ndir value.")
			sys.exit()
