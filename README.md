# AkaiKKRPythonUtil
Python utilities for AkaiKKR

# License

 Copyright (c) 2021-2023 AkaiKKRteam.
 Distributed under the terms of the Apache License, Version 2.0.



# installation
{PREFIX} is the directory you installed AkaiKKRPythonUtil.

install PyAkaiKKR by
```
cd {PREFIX}/library/PyAkaiKKR
pip install .
```

If you use {PREFIX}/tests, also install AkaiKKRTestScript by
```
cd {PREFIX}/library/AkaiKKRTestScript
pip install .
```

## for the CPA2021V01 user

{AKAIKKRCPA2021V01} is the directory where you installed CPA2021V01.
{AKAIKKRCPA2021V01} can be an absolute or a relative path.

1. Please apply patch to akaikkr_cpa2021v01.
```
$ cd {AKAIKKRCPA2021V01}
$ patch -p1 < {PREFIX}/tests/akaikkr_cpa2021v01/akaikkr_cpa2021v01.patch
```

The patched cpa2021v01 has 
- 'geom' mode
- Site labels exteded to 40 characters.
- additional information on the "dispersion" files.

2. make specx and fmg in CPA2021V01.
```
$ make specx fmg
```

# test script

- It tests AkaiKKR package.
- It generates AkaiKKR input file from cif files.

 
For example, if you use CPA2021V01, (Note that the directory name specified by {PREFIX} must be .../akaikkr_cpa2021v01.)
```
$ cd {PREFIX}
$ python testrun.py
```

The following display appears at the end of the execution.
```
SHORT SUMMARY
             dos fsm go gofmg spc31
AlMnFeCo_bcc   O   O  O     O     -
Co             O   O  O           -
Co2MnSi        O   O  O           -
Cu             O      O           -
Fe             O   O  O           -
FeB195         O      O           -
FeRh05Pt05     O   O  O           -
Fe_lmd         O      O           -
GaAs           O      O           -
Ni             O   O  O           -
NiFe           O   O  O           -
SmCo5_noc             O            
SmCo5_oc       O   O  O           -
O: passed. X: failed, -: no reference
```

As an example, the block spectra $A(w,k)$ of NiFe (FCC $\mathrm{Ni}_{0.9}\mathrm{Fe}_{0.1}$) is generated under the NiFe directory as follows.
![](https://github.com/nim-hrkn/AkaiKKRPythonUtil/blob/cpa2021v01_supported/fig/NiFe_Awk_all.png?raw=true)

# BUG
- TEST FAILED is always shown at the end of testrun.py.
- Awk\_both.png is generated at the top directory of testrun.py.
