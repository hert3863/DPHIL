#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 09:24:23 2018

@author: bakerh
"""

import numpy as np


def cut(Ni=64, T=5, iterate=10):
    dist = np.zeros(iterate)
    # perform cut iteration 'iterate' number of times
    for i in range(iterate):
        # initialise length
        N = Ni
        # iterate T times
        for j in range(T):
            # number of interior coordinates
            intN = N-1
            # if N is 2 or less then can't cut twice uniquely
            if N > 2:
                # random cut locations
                c1 = np.random.randint(1, intN+1)
                c2 = np.random.randint(1, intN+1)
                # ensure cut locations are unique
                while c1 == c2:
                    c2 = np.random.randint(1, intN+1)
                # ensure cut loc 1 is lower than loc 2
                if c1 > c2:
                    hold = c1
                    c1 = c2
                    c2 = hold
                # compute lengths of 3 new ropes
                l1 = c1
                l2 = c2-c1
                l3 = N-c2
                # select longest rope
                N = np.max([l1, l2, l3])
            else:
                break
        dist[i] = N
    return dist


# compute distributions for differing initial conditions
dist64 = cut(Ni=64, T=5, iterate=10000000)
dist1024 = cut(Ni=1024, T=10, iterate=10000000)

a1 = np.mean(dist64)
a2 = np.std(dist64)
a3 = np.mean(dist1024)
a4 = np.std(dist1024)

dist4 = dist64[dist64 > 4]
dist8 = dist4[dist4 > 8]

a5 = len(dist8)/len(dist4)

dist6 = dist1024[dist1024 > 6]
dist12 = dist6[dist6 > 12]

a6 = len(dist12)/len(dist6)






















