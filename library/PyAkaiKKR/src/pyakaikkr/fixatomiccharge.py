from pymatgen.core import Structure


def fix_charged_structure(structure: Structure):
    """ fix atomic chage.

    pymatgen can't handle atomic charges, such as Fe2.33+, As3-, K+,
    The atomic charges will be removed in this function.

    Args:
        structure (pymatgen.core.Structure): periodic structure.

    Returns:
        pymatgen.core.Structure: atomic-charge-fixed structure.
    """
    import re
    replace_site_list = []
    for _i, site in enumerate(structure.sites):
        species = site.species
        replace_dic = {}
        found = False
        for _elmname in species.as_dict().keys():
            _replaced_elmname = re.sub(r"[0-9+-].*$", "", _elmname)
            replace_dic.update({_elmname: _replaced_elmname})
            if _elmname != _replaced_elmname:
                found = True
        if found:
            species = species.replace(replace_dic)
            replace_site_list.append([_i, species])
    for _i, species in replace_site_list:
        structure.replace(_i, species)
    return structure
