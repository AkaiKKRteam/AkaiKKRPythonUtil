# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.


import numpy as np
from typing import List
import os
import subprocess
from math import atan, cos, sin, pi, sqrt, fabs

if False:
    def get_prmvec(brav: str, coa: float, boa: float, alpha: float, beta: float, gamma: float):

        prmrconv = PrmrConverter()
        prmr = prmrconv.get(brav)
        if brav == "hcp":
            """
          boa=1e0
          if(abs(coa) < 1e-6) coa=sqrt(8d0/3d0)
          alpha=90d0
          beta=90d0
          gamma=120d0
          go to 170
        """
            boa = 1.0
            if abs(coa) < 1e-6:
                cos = math.sqrt(8.0/3.0)
                alpha = 90.0
                beta = 90.0
                gamma = 12.0
        else:
            raise ValueError("bravtyp={} not supported".format(brav))

        """ 
        do 210 i=1,3
        r(1,i)=prmr(1,i,ibrav)
        r(2,i)=prmr(2,i,ibrav)*boa
    210 r(3,i)=prmr(3,i,ibrav)*coa
        """
        prmr[:, 1] = prmr[:, 1] * boa
        prmr[:, 2] = prmr[:, 2] * coa

        return prmr


class PrmrConverter:
    """mimic akaikkr/source/prmvec.f and ibrav.f and vrotat.f
    """

    def __init__(self):
        c0 = 0.0
        c1 = 0.5
        c2 = -0.5
        c3 = 1.0
        c4 = -1
        c5 = 8.660254037844386e-1
        c5 = cos(pi*1.0/6.0)
        c6 = -c5

        _prmr0 = np.zeros(9).tolist()
        _prmr = [c0, c1, c1, c1, c0, c1, c1, c1, c0,  c2, c1, c1, c1, c2, c1, c1, c1, c2,
                 c1, c6, c0, c1, c5, c0, c0, c0, c3,  c3, c0, c0, c0, c3, c0, c0, c0, c3,
                 c2, c1, c1, c1, c2, c1, c1, c1, c2,  c3, c0, c0, c0, c3, c0, c0, c0, c3,
                 c0, c1, c1, c1, c0, c1, c1, c1, c0,  c2, c1, c1, c1, c2, c1, c1, c1, c2,
                 c1, c2, c0, c1, c1, c0, c0, c0, c3,  c3, c0, c0, c0, c3, c0, c0, c0, c3,
                 c1, c2, c0, c1, c1, c0, c0, c0, c3,  c3, c0, c0, c0, c3, c0, c0, c0, c3,
                 c3, c0, c0, c0, c3, c0, c0, c0, c3,  c0, c0, c0, c0, c0, c0, c0, c0, c0,
                 c0, c1, c1, c1, c0, c1, c1, c1, c0]

        _prmr0.extend(_prmr)

        self.prmr = np.array(_prmr0).reshape(16, 3, 3)

    def ibrava(self, brvtyp: str) -> int:
        """
         from source/ibrava.f

        data brav/'fcc','bcc','hcp','sc','bct','st','fco','bco',
     &          'bso','so','bsm','sm','trc','rhb','fct','trg','hex',
     &          'aux','prv'/,
     &     ibrav/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,14,3,16,16/
        """
        bravlist = ['fcc', 'bcc', 'hcp', 'sc', 'bct', 'st', 'fco', 'bco',
                    'bso', 'so', 'bsm', 'sm', 'trc', 'rhb', 'fct', 'trg', 'hex',
                    'aux', 'prv']
        ibrav = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                 11, 12, 13, 14, 15, 14, 3, 16, 16]

        i = bravlist.index(brvtyp)
        return ibrav[i]

    def apply170(self, ibrav: int, coa: float, boa: float,
                 angl: List[float] = None) -> np.ndarray:
        """mimic label 170

            do 210 i=1,3
            r(1,i)=prmr(1,i,ibrav)
            r(2,i)=prmr(2,i,ibrav)*boa
        210 r(3,i)=prmr(3,i,ibrav)*coa
        """
        r = self.prmr[ibrav]
        r[:, 1] = r[:, 1]*boa
        r[:, 2] = r[:, 2]*coa

        if angl is None:
            angl = [0.0, 0.0, 0.0]
        return self.apply180(r, angl)

    def apply180(self, r: np.ndarray, angl: List[float]) -> np.ndarray:
        """ mimic label 180 and vrotat.f
        c--     (v1(1),v1(2),v1(3)) -> (t,y,v1(3))
        c=cos(alpha)
        s=sin(alpha)
        t=c*v1(1)-s*v1(2)
        y=s*v1(1)+c*v1(2)
  c---    (t,y,v1(3)) -> (x,y,v2(3))
        c=cos(beta)
        s=sin(beta)
        v2(3)=c*v1(3)-s*t
        x=s*v1(3)+c*t
  c---     (x,y,v2(3)) -> (v2(1),v2(2),v2(3))
        c=cos(gamma)
        s=sin(gamma)
        v2(1)=c*x-s*y
        v2(2)=s*x+c*y
      """
        rad = pi/180.0
        # order is 2,1,0
        alpha = angl[2]*rad
        beta = angl[1]*rad
        gamma = angl[0]*rad

        v2all = np.zeros_like(r)

        for i, _v1 in enumerate(r):
            v1 = np.zeros((4))
            v1[1:] = _v1

            v2 = np.zeros((4))

            c = cos(alpha)
            s = sin(alpha)
            t = c*v1[1]-s*v1[2]
            y = s*v1[1]+c*v1[2]
            c = cos(beta)
            s = sin(beta)
            v2[3] = c*v1[3]-s*t
            x = s*v1[3]+c*t
            c = cos(gamma)
            s = sin(gamma)
            v2[1] = c*x-s*y
            v2[2] = s*x+c*y

            v2all[i, 0:3] = v2[1:4]

        return v2all

    def apply(self, brvtyp: str, coa: float, boa: float,
              alpha: float, beta: float, gamma: float,
              r: np.ndarray = None,
              angl: List[float] = None) -> np.ndarray:
        """calculate matrix r.

        from source/prmvec.f

        Args:
            brvtyp (str): brvtype
            coa (float): c/a
            boa (float): b/a
            alpha (float): alpha angle
            beta (float): beta angle
            gamma (float): gamma angle
            r (np.ndarray, optional): in the case of prv or aux. Defaults to None.
            angl (List[float, float, float], optional): distortion angle. Defaults to None.

        Returns:
            np.ndarray: primitive vector
        """

        if angl is None:
            angl = np.array([0.0, 0.0, 0.0])

        ibrav = self.ibrava(brvtyp)
        if r is not None:
            r = np.array(r)
#      go to (10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,160),
#     &       ibrav
        if ibrav == 1:
            #     --- fcc case ---
            boa = 1e0
            coa = 1e0
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 2:
            #     --- bcc case ---
            boa = 1e0
            coa = 1e0
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 3:
            #     --- hcp case ---
            boa = 1e0
            if(abs(coa) < 1e-6):
                coa = sqrt(8e0/3e0)
            alpha = 90e0
            beta = 90e0
            gamma = 120e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 4:
            #     --- sc case ---
            boa = 1e0
            coa = 1e0
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 5:
            #     --- bct (fct) case ---
            boa = 1e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 6:
            #     --- st case ---
            boa = 1e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 7:
            #     --- fco case ---
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 8:
            #     --- bco case ---
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 9:
            #     --- bso ---
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 10:
            #     --- so case ---
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            # go to 170
            r = self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 11:
            #     --- bsm case ---
            #  110 continue
            alpha = 90e0
            gamma = 90e0
            if(abs(beta) < 1e-6):
                beta = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            r = [[5e-1, -5e-1*boa, 0e0],
                 [5e-1, 5e-1*boa, 0e0],
                 [cos(pi/180e0*beta)*coa, 0e0, sqrt(1e0-(cos(pi/180e0*beta))**2)*coa]]
            # go to 180
            return self.apply180(r, angl)
        elif ibrav == 12:
            #     --- sm case ---
            #  120 continue
            alpha = 90e0
            gamma = 90e0
            if(abs(beta) < 1e-6):
                beta = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            r = [[1e0, 0e0, 0e0],
                 [0e0, boa, 0e0],
                 [cos(pi/180e0*beta)*coa, 0e0, sqrt(1e0-(cos(pi/180e0*beta))**2)*coa]]
            # go to 180
            return self.apply180(r, angl)
        elif ibrav == 13:
            #     --- trc case ---
            #  130 continue
            if(abs(alpha) < 1e-6):
                alpha = 90e0
            if(abs(beta) < 1e-6):
                beta = 90e0
            if(abs(gamma) < 1e-6):
                gamma = 90e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            if(abs(boa) < 1e-6):
                boa = 1e0
            cw = 1e0-(cos(pi/180e0*alpha))**2-(cos(pi/180e0*beta))**2 -\
                (cos(pi/180e0*gamma))**2+2*cos(pi/180e0*alpha) *\
                cos(pi/180e0*beta)*cos(pi/180e0*gamma)
            r = [[1e0, 0e0, 0e0],
                 [cos(pi/180e0*gamma)*boa, sin(pi/180e0*gamma)*boa, 0e0],
                 [cos(pi/180e0*beta)*coa,
                 (cos(pi/180e0*alpha)-cos(pi/180e0*beta) *
                  cos(pi/180e0*gamma))/sqrt(1e0-(cos(pi/180e0*gamma))**2)*coa,
                  sqrt(cw)/sqrt(1e0-(cos(pi/180e0*gamma))**2)*coa]]
            # go to 180
            return self.apply180(r, angl)
        elif ibrav == 14:
            #     --- rhb (trg) case ---
            #     --- obverse setting is used
            #         Note: in the previous version, the reverse setting is used.
            boa = 1e0
            coa = 1e0
            if(abs(alpha) < 1e-6):
                alpha = 60e0
            beta = alpha
            gamma = alpha
            theta = 2e0*pi*alpha/360e0
            z = sqrt((5e-1+cos(theta))/(1e0-cos(theta)))
            snrm = sqrt(1e0+z**2)
            r = [[(sqrt(3e0)/2e0)/snrm, -(1e0/2e0)/snrm, z/snrm],
                 [0e0, 1e0/snrm, z/snrm],
                 [-(sqrt(3e0)/2e0)/snrm, -(1e0/2e0)/snrm, z/snrm]]
            # go to 180
            return self.apply180(r, angl)
        elif ibrav == 15:
            #     --- fct (bct) case ---
            boa = 1e0
            if(abs(coa) < 1e-6):
                coa = 1e0
            alpha = 90e0
            beta = 90e0
            gamma = 90e0
            # go to 170
            return self.apply170(ibrav, coa, boa, angl)
        elif ibrav == 16 or ibrav == 17:
            #  aux or prv
            """
            do 190 i=1,3
        190 rlng(i)=sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2)
            """
            rlng = []
            for r1 in r:
                rlng.append(np.sqrt(np.sum(r1**2)))

            rlng = np.array(rlng)
            """
          do 200 i=1,3
          j=mod(i,3)+1
      200 rprd(i)=(r(1,i)*r(1,j)+r(2,i)*r(2,j)+r(3,i)*r(3,j))
         &        /rlng(i)/rlng(j)
            """
            rprd = []
            for i in range(3):
                j = (i+1) % 3
                # (ri* rj)/|ri|{rj| = cos
                rprd.append((r[i, 0]*r[j, 0] + r[i, 1]*r[j, 1] +
                            r[i, 2]*r[j, 2])/rlng[i]/rlng[j])
            rprd = np.array(rprd)

            zero = 1e-6
            """
            if(    abs(rlng(1)-rlng(2)) < zero
          & .and. abs(rlng(2)-rlng(3)) < zero
      #     --- a=b=c
          & .and. abs(rprd(1)-rprd(2)) < zero
          & .and. abs(rprd(2)-rprd(3)) < zero
      #     --- alpha=beta=gamma
          & .and. abs(r(3,1)-r(3,2)) < zero
          & .and. abs(r(3,2)-r(3,3)) < zero) then
      #     --- z-component of a, b, c is zero
      #     --- this case should be trigonal ---
            """

            if fabs(rlng[0]-rlng[1]) < zero and \
                    fabs(rlng[1]-rlng[2]) < zero and \
                    fabs(rprd[0]-rprd[1]) < zero and \
                    fabs(rprd[1]-rprd[2]) < zero and  \
                    fabs(r[0, 2]-r[1, 2]) < zero and fabs(r[1, 2]-r[2, 2]) < zero:
                ibrav = 17
                coa = 1e0
                boa = 1e0
                # go to 180
                return self.apply180(r, angl)
                # endif
            """"
            if(    abs(rlng(1)-rlng(2)) < zero
      #     --- a=b
          &     .and. (abs(rprd(1)-5e-1) < zero
      #     --- gamma=60
          &     .or. abs(rprd(1)+5e-1) <.zero)
      #     --- gamma=120
          &     .and. abs(rprd(2)) < zero
          &     .and. abs(rprd(3)) < zero) then
      #     --- alpha=beta=90
      #     --- this case should be hcp ---
          """
            if fabs(rlng[0]-rlng[1]) < zero and \
                    fabs(rprd[0]-5e-1) < zero or  \
                    fabs(rprd[0]+5e-1) < zero and  \
                    fabs(rprd[1]) < zero and fabs(rprd[2]) < zero:
                #     write(*,'(1x,a)')'this case is hexagonal'
                ibrav = 17
                coa = 2e0*r[2, 2]
                boa = 1e0
                # go to 180
                return self.apply180(r, angl)
                # endif
    #     --- at the moment not fixed but likeky neither hex nor trg---
        ibrav = 16
        # go to 180
        return self.apply180(r, angl)

        raise ValueError("unknown brvtype={}".format(brvtype))


if __name__ == "__main__":
    class PrmrRun:
        """execute a simple program to check prmvec
        """

        def __init__(self, directory=".", inputfile="input.txt",
                     outputfile="output.txt", exe="./prmvectest"):
            self.directory = directory
            self.inputfile = inputfile
            self.outputfile = outputfile
            self.exe = exe

        def set_param(self, brvtyp: str, coa: float, boa: float,
                      alpha: float, beta: float, gamma: float, r: np.ndarray, angl: np.ndarray):
            """write parameters to the file given by self.inputfile

            Args:
                brvtyp (str): [description]
                coa (float): [description]
                boa (float): [description]
                alpha (float): [description]
                beta (float): [description]
                gamma (float): [description]
                r (np.ndarray): [description]
                angl (np.ndarray): 
            """
            if r is None:
                r = np.array([0.0 for i in range(9)])
            else:
                r = np.array(r)
                r = r.ravel()
            if angl is None:
                angl = np.array([0.0 for i in range(3)])
            else:
                angl = np.array(angl)

            inputfile = os.path.join(self.directory, self.inputfile)
            with open(inputfile, "w") as f:
                p = [brvtyp, coa, boa, alpha, beta, gamma]
                p.extend(r.tolist())
                p.extend(angl.tolist())
                p = list(map(str, p))
                f.write("\n".join(p))

        def get_result(self) -> np.ndarray:
            """run program and get result

            Returns:
                np.ndarray: r matrix
            """
            cmd = "(cd {}; {} < {} > {})".format(self.directory, self.exe,
                                                 self.inputfile,  self.outputfile)
            subprocess.call(cmd, shell=True)
            outputfile = os.path.join(self.directory, self.outputfile)
            with open(outputfile) as f:
                data = f.read()

            s = data.split()
            s = list(map(float, s))
            r = np.array(s).reshape(3, 3)
            return r

    def compare(p):
        """compare python code and fortran code

        Args:
            p (dict): kkr parameter
        """
        brvtyp = p["brvtyp"]
        coa = p["coa"]
        boa = p["boa"]
        alpha = p["alpha"]
        beta = p["beta"]
        gamma = p["gamma"]
        if "r" in p:
            r = np.array(p["r"])
        else:
            r = None
        if "angl" in p:
            angl = np.array(p["angl"])
        else:
            angl = None

        # prmr = get_prmvec("hcp", coa, boa, alpha, beta, gamma)
        prmrconv = PrmrConverter()
        prmr = prmrconv.apply(brvtyp, coa, boa, alpha, beta, gamma, r, angl)
        print("python", prmr)
        np.set_printoptions(precision=15)

        prmrrun = PrmrRun()
        prmrrun.set_param(brvtyp, coa, boa, alpha, beta, gamma, r, angl)
        runr = prmrrun.get_result()
        print("runr", runr)

        eps = 1e-10
        if np.all(np.abs(prmr.ravel() - runr.ravel()) < eps):
            print(brvtyp, "same.")
        else:
            print(brvtyp, "different.")
            print("python", prmr)
            print("runr", runr)
            raise ValueError("comparison failed")

    def main():
        coa = 0.0
        boa = 0.0
        alpha = 0.0
        beta = 0.0
        gamma = 0.0
        if True:
            p = {"brvtyp": "hcp",
                 "coa": coa,
                 "boa": boa,
                 "alpha": alpha,
                 "beta": beta,
                 "gamma": gamma}
            compare(p)

            p = {"brvtyp": "fcc",
                 "coa": coa,
                 "boa": boa,
                 "alpha": alpha,
                 "beta": beta,
                 "gamma": gamma}
            compare(p)

            p = {"brvtyp": "bcc",
                 "coa": coa,
                 "boa": boa,
                 "alpha": alpha,
                 "beta": beta,
                 "gamma": gamma}
            compare(p)

            p = {"brvtyp": "sc",
                 "coa": coa,
                 "boa": boa,
                 "alpha": alpha,
                 "beta": beta,
                 "gamma": gamma}
            compare(p)

        p = {"brvtyp": "bct",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "st",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "fco",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "bco",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "bso",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "so",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "bsm",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "sm",
             "coa": coa,
             "boa": boa,

             "alpha": alpha,
             "beta": beta,
             "gamma": 30.0}
        compare(p)

        p = {"brvtyp": "trc",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "rhb",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": 30.0}
        compare(p)

        p = {"brvtyp": "fct",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "trg",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        p = {"brvtyp": "hex",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma}
        compare(p)

        # spacial case
        r = [[0.642787609686539, -0.37111359948428, 0.67014833067858],
             [0.,   0.742227198968559, 0.67014833067858],
             [-0.642787609686539, -0.37111359948428, 0.67014833067858]]

        p = {"brvtyp": "prv",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma,
             "r": r}

        compare(p)

        r = [[0.5,    -0.866025403784439, 0.],
             [0.5,  0.866025403784439, 0.],
             [0.,  0.,  1.621513944223108]]

        p = {"brvtyp": "aux",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma,
             "r": r}
        compare(p)

        angl = [10.0, 5.0, 6.0]
        p = {"brvtyp": "aux",
             "coa": coa,
             "boa": boa,
             "alpha": alpha,
             "beta": beta,
             "gamma": gamma,
             "r": r,
             "angl": angl}
        compare(p)

    main()
