# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 16:31:53 2017

@author: bakerh
"""

import numpy as np


def jetindicesmean(u, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    iGCM output

    Parameters
    ----------
    u: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    indices: float
        max jet speed and position of maximum
        in winter/summer, speed/lat format
    """
    u850 = np.copy(u)
    indices = np.zeros((np.ma.size(u850, axis=0), 2))
    u850 = np.mean(u850[:, 26:94, 400:], axis=2)  # 147:334
    latna = lat[26:94]
    for t in range(np.ma.size(u850, axis=0)):
        indices[t, 0] = np.amax(u850[t, :])
        indices[t, 1] = latna[np.argmax(u850[t, :])]
    return indices


def lanczos(unfiltered):
    hold = []
    W = [0.446640113241E-08,
         -0.274115602902E-02,
         -0.461534097037E-02,
         -0.202356429576E-07,
         0.964742314469E-02,
         0.129707987317E-01,
         0.374235240357E-07,
         -0.219366594221E-01,
         -0.281014728087E-01,
         -0.528782831460E-07,
         0.466935604307E-01,
         0.620463289111E-01,
         0.635961703709E-07,
         -0.134316831363E+00,
         -0.273896602685E+00,
         0.666666600000E+00,
         -0.273896602685E+00,
         -0.134316831363E+00,
         0.635961703709E-07,
         0.620463289111E-01,
         0.466935604307E-01,
         -0.528782831460E-07,
         -0.281014728087E-01,
         -0.219366594221E-01,
         0.374235240357E-07,
         0.129707987317E-01,
         0.964742314469E-02,
         -0.202356429576E-07,
         -0.461534097037E-02,
         -0.274115602902E-02,
         0.446640113241E-08]

    for n in range(13956-30):
        sums = 0
        for m in range(n, n+31):
            a = unfiltered[m] * W[m-n]
            sums = sums + a
        hold = np.append(hold, sums)
    filtered = hold
    return filtered
