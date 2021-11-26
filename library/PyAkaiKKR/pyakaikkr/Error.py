# coding: utf-8
# Copyright (c) 2021 AkaiKKRteam.
# Distributed under the terms of the Apache License, Version 2.0.



class CIF2KKRGetStructureError(Exception):
    """failed to get_structure(primitive=primitive)"""
    pass


class CIF2KKRGetConventionalStandardStructureError(Exception):
    """failed to get_conventional_standard_structure"""
    pass


class CIF2KKRSpgDifferentError(Exception):
    """spg number in the cif file != spg number from spglib"""
    pass


class CIF2KKRCellShapeError(Exception):
    """invalid cell shape. inconsistent with the definition of AkaiKKR"""
    pass


class CIF2KKRUnknownElementError(Exception):
    """unknown element. Proabably Element name isn't defined correctly."""
    pass


class CIF2KKRNoStructureError(Exception):
    """no strcuture is read."""
    pass


class CIF2KKRTooManyTypesError(Exception):
    """to many types are defines"""
    pass


class CIF2KKRNsiteInconsistentError(Exception):
    """number of sites in cif primitive isn't equal to the number of KKR sites"""
    pass


class CIF2KKRTooManyElementError(Exception):
    """number of sites in cif primitive isn't equal to the number of KKR sites"""
    pass


class KKRFailedExecutionError(Exception):
    """failed to execute specx"""
    pass


class KKRValueAquisitionError(Exception):
    """failed to execute specx"""
    pass
