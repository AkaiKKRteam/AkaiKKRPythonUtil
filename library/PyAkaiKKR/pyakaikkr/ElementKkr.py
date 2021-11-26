# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.

from copy import deepcopy


class ElementKKR:
    """define default values of Element."""
    class _Element:
        def __init__(self, number, rmt, field, lmax):
            self.number = number
            self.rmt = rmt
            self.field = field
            self.lmax = lmax

        def __str__(self):
            return "Z={} rmt={} field={} lmax={}".format(self.number,
                                                         self.rmt,
                                                         self.field,
                                                         self.lmax)

    _dict = {  # Z, rmt, field, lmax
        # 'Vc': _Element(0, 0.0, 0.0, 2),
        'H': _Element(1, 0.0, 0.0, 2),
        'He': _Element(2, 0.0, 0.0, 2),
        'Li': _Element(3, 0.0, 0.0, 2),
        'Be': _Element(4, 0.0, 0.0, 2),
        'B': _Element(5, 0.0, 0.0, 2),
        'C': _Element(6, 0.0, 0.0, 2),
        'N': _Element(7, 0.0, 0.0, 2),
        'O': _Element(8, 0.0, 0.0, 2),
        'F': _Element(9, 0.0, 0.0, 2),
        'Ne': _Element(10, 0.0, 0.0, 2),
        'Na': _Element(11, 0.0, 0.0, 2),
        'Mg': _Element(12, 0.0, 0.0, 2),
        'Al': _Element(13, 0.0, 0.0, 2),
        'Si': _Element(14, 0.0, 0.0, 2),
        'P': _Element(15, 0.0, 0.0, 2),
        'S': _Element(16, 0.0, 0.0, 2),
        'Cl': _Element(17, 0.0, 0.0, 2),
        'Ar': _Element(18, 0.0, 0.0, 2),
        'K': _Element(19, 0.0, 0.0, 2),
        'Ca': _Element(20, 0.0, 0.0, 2),
        'Sc': _Element(21, 0.0, 0.0, 2),
        'Ti': _Element(22, 0.0, 0.0, 2),
        'V': _Element(23, 0.0, 0.0, 2),
        'Cr': _Element(24, 0.0, 0.0, 2),
        'Mn': _Element(25, 0.0, 0.0, 2),
        'Fe': _Element(26, 0.0, 0.0, 2),
        'Co': _Element(27, 0.0, 0.0, 2),
        'Ni': _Element(28, 0.0, 0.0, 2),
        'Cu': _Element(29, 0.0, 0.0, 2),
        'Zn': _Element(30, 0.0, 0.0, 2),
        'Ga': _Element(31, 0.0, 0.0, 2),
        'Ge': _Element(32, 0.0, 0.0, 2),
        'As': _Element(33, 0.0, 0.0, 2),
        'Se': _Element(34, 0.0, 0.0, 2),
        'Br': _Element(35, 0.0, 0.0, 2),
        'Kr': _Element(36, 0.0, 0.0, 2),
        'Rb': _Element(37, 0.0, 0.0, 2),
        'Sr': _Element(38, 0.0, 0.0, 2),
        'Y': _Element(39, 0.0, 0.0, 2),
        'Zr': _Element(40, 0.0, 0.0, 2),
        'Nb': _Element(41, 0.0, 0.0, 2),
        'Mo': _Element(42, 0.0, 0.0, 2),
        'Tc': _Element(43, 0.0, 0.0, 2),
        'Ru': _Element(44, 0.0, 0.0, 2),
        'Rh': _Element(45, 0.0, 0.0, 2),
        'Pd': _Element(46, 0.0, 0.0, 2),
        'Ag': _Element(47, 0.0, 0.0, 2),
        'Cd': _Element(48, 0.0, 0.0, 2),
        'In': _Element(49, 0.0, 0.0, 2),
        'Sn': _Element(50, 0.0, 0.0, 2),
        'Sb': _Element(51, 0.0, 0.0, 2),
        'Te': _Element(52, 0.0, 0.0, 2),
        'I': _Element(53, 0.0, 0.0, 2),
        'Xe': _Element(54, 0.0, 0.0, 2),
        'Cs': _Element(55, 0.0, 0.0, 2),
        'Ba': _Element(56, 0.0, 0.0, 2),
        'La': _Element(57, 0.0, 0.0, 2),
        'Ce': _Element(58, 0.0, 0.0, 2),
        'Pr': _Element(59, 0.0, 0.0, 2),
        'Nd': _Element(60, 0.0, 0.0, 2),
        'Pm': _Element(61, 0.0, 0.0, 2),
        'Sm': _Element(62, 0.0, 0.0, 2),
        'Eu': _Element(63, 0.0, 0.0, 2),
        'Gd': _Element(64, 0.0, 0.0, 2),
        'Tb': _Element(65, 0.0, 0.0, 2),
        'Dy': _Element(66, 0.0, 0.0, 2),
        'Ho': _Element(67, 0.0, 0.0, 2),
        'Er': _Element(68, 0.0, 0.0, 2),
        'Tm': _Element(69, 0.0, 0.0, 2),
        'Yb': _Element(70, 0.0, 0.0, 2),
        'Lu': _Element(71, 0.0, 0.0, 2),
        'Hf': _Element(72, 0.0, 0.0, 2),
        'Ta': _Element(73, 0.0, 0.0, 2),
        'W': _Element(74, 0.0, 0.0, 2),
        'Re': _Element(75, 0.0, 0.0, 2),
        'Os': _Element(76, 0.0, 0.0, 2),
        'Ir': _Element(77, 0.0, 0.0, 2),
        'Pt': _Element(78, 0.0, 0.0, 2),
        'Au': _Element(79, 0.0, 0.0, 2),
        'Hg': _Element(80, 0.0, 0.0, 2),
        'Tl': _Element(81, 0.0, 0.0, 2),
        'Pb': _Element(82, 0.0, 0.0, 2),
        'Bi': _Element(83, 0.0, 0.0, 2),
        'Po': _Element(84, 0.0, 0.0, 2),
        'At': _Element(85, 0.0, 0.0, 2),
        'Rn': _Element(86, 0.0, 0.0, 2),
        'Fr': _Element(87, 0.0, 0.0, 2),
        'Ra': _Element(88, 0.0, 0.0, 2),
        'Ac': _Element(89, 0.0, 0.0, 2),
        'Th': _Element(90, 0.0, 0.0, 2),
        'Pa': _Element(91, 0.0, 0.0, 2),
        'U': _Element(92, 0.0, 0.0, 2),
        'Np': _Element(93, 0.0, 0.0, 2),
        'Pu': _Element(94, 0.0, 0.0, 2),
        'Am': _Element(95, 0.0, 0.0, 2),
        'Cm': _Element(96, 0.0, 0.0, 2),
        'Bk': _Element(97, 0.0, 0.0, 2),
        'Cf': _Element(98, 0.0, 0.0, 2),
        'Es': _Element(99, 0.0, 0.0, 2),
        'Fm': _Element(100, 0.0, 0.0, 2),
        'Md': _Element(101, 0.0, 0.0, 2),
        'No': _Element(102, 0.0, 0.0, 2),
        'Lr': _Element(103, 0.0, 0.0, 2),
        'Rf': _Element(104, 0.0, 0.0, 2),
        'Db': _Element(105, 0.0, 0.0, 2),
        'Sg': _Element(106, 0.0, 0.0, 2),
        'Bh': _Element(107, 0.0, 0.0, 2),
        'Hs': _Element(108, 0.0, 0.0, 2),
        'Mt': _Element(109, 0.0, 0.0, 2),
        'Ds': _Element(110, 0.0, 0.0, 2),
        'Rg': _Element(111, 0.0, 0.0, 2),
        'Cn': _Element(112, 0.0, 0.0, 2),
        'Nh': _Element(113, 0.0, 0.0, 2),
        'Fl': _Element(114, 0.0, 0.0, 2),
        'Mc': _Element(115, 0.0, 0.0, 2),
        'Lv': _Element(116, 0.0, 0.0, 2),
        'Ts': _Element(117, 0.0, 0.0, 2),
        'Og': _Element(118, 0.0, 0.0, 2),
    }

    def __init__(self, Vc):
        """constructure

        Args:
            Vc (str, optional): element treated as Z=0. Defaults to None
        """
        self.dict = deepcopy(self._dict)
        if Vc is not None:
            elm = self.dict[Vc]
            self.dict[Vc] = self._Element(
                0, elm.rmt, elm.field, elm.lmax)
            print("debug, _ElementKKR", Vc, ":", self.dict[Vc])

    def getAtomicNumber(self, element):
        if element in self.dict:
            return self.dict[element].number
        else:
            return 0

    def getAtomicRMT(self, element):
        if element in self.dict:
            return self.dict[element].rmt
        else:
            return 0.0

    def getAtomicField(self, element):
        if element in self.dict:
            return self.dict[element].field
        else:
            return 0.0

    def getAtomicLMax(self, element):
        if element in self.dict:
            return self.dict[element].lmax
        else:
            return 0
