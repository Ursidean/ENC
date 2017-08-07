"""
This module calculate the area-weighted absolute average clumpiness error on
a per class basis across the actively allocated land-use classes.
"""

from clumpy_module import clumpiness_index


def AACE(amap, smap, mask, pas, act, luc, luc_size):
    # Calculate the clumpiness for the actual map.
    amap_clu = clumpiness_index(amap, mask, luc)
    # Calculate the clumpiness for the simulated map.
    smap_clu = clumpiness_index(smap, mask, luc)
    # Calculate the absolute error across the active classes.
    abs_class_clu_error = [0]*act
    for i in range(0, act):
        abs_class_clu_error[i] = abs(amap_clu[i + pas] - smap_clu[i + pas])
    # Now calculate the weighted sum of all
    AW_ACE = 0
    for i in range(0, act):
        AW_ACE = AW_ACE + (luc_size[i + pas]*abs_class_clu_error[i])
    AW_ACE = AW_ACE/sum(luc_size)
    return AW_ACE
