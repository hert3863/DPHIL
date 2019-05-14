#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 13:57:05 2018

@author: bakerh
"""

import numpy as np


def eof_analysis(control, lat, lon, response, region=[-90, 90, 0, 360],
                 neof=20):
    '''
    Amalgamates all EOF functions into one process

    Parameters
    ----------
    control: array
        data to determine EOFs for. time x nlat x nlon
    lat: array
        full data latitude
    lon: array
        full data longitude
    response: array
        response data to project onto eofs. n_forceruns x nlat x nlon
    region: gridpoints
        region to compute eofs over
    nth: int
        number of eofs to compute

    Returns
    -------
    eofs: array
        First n EOFs. neof x nlat x nlon
    var: array
        amount of variance explained by each EOF. neof x 3 array:
        var[:, 0] = the 1st neof eigenvalues of the covariance matrix
        var[:, 1] = the 1st neof explained variances (in %).
        var[:, 2] = the 1st neof 1sd errors (in %).
    proj: array
        projection of response field onto the first n EOFs
    '''
    # compute eofs

    eofs_un, pcs, var = eof_svd(control, lat, lon, neof, region)
    # regress control onto pcs to get sensible units
    eofs = np.zeros((neof, len(lat), len(lon)))
    for i in range(neof):
        eofs[i] = eof_regress(pcs, i, control.copy())
    # project responses onto eof
    proj = np.zeros((neof, len(response)))

    for i in range(neof):
        proj[i, :] = eof_response(response, eofs[i], lat, lon, region)
    return eofs, var, proj


def eof_svd(control, lati, long, neof, region=[-90, 90, 0, 360]):
    '''
    EOFs of data

    Parameters
    ----------
    control: array
        data to determine EOFs for. time x nlat x nlon


    Returns
    -------
    eofs: array
        First n EOFs. neof x nsigma x nlat
    var: array
        amount of variance explained by each EOF. neof x 3 array:
        var[:, 0] = the 1st neof eigenvalues of the covariance matrix
        var[:, 1] = the 1st neof explained variances (in %).
        var[:, 2] = the 1st neof errors (in %).
    pcs: array
        principal components. time x neof
    '''
    from numpy import linalg
    field = control.copy()[:, (region[0] <= lati) & (lati <= region[1]), :]
    field = field[:, :, (region[2] <= long) & (long <= region[3])]
    lat = lati[(region[0] <= lati) & (lati <= region[1])]
    lon = long[(region[2] <= long) & (long <= region[3])]

    dim = field.shape
    field = eofweight(field, lat, lon)
    aver = np.mean(field, axis=0)
    for i in range(dim[0]):
        field[i, :] = field[i, :] - aver
    field = np.reshape(field, (dim[0], dim[1] * dim[2]), order='F')

    if dim[0] > dim[1] * dim[2]:
        field = np.transpose(field)

    u, s, v = linalg.svd(field)
    v = np.transpose(v)
    s = s * s
    eigs = np.copy(s)
    s = s / np.sum(s)
    if dim[0] < dim[1] * dim[2]:
        u, v = v, u

    eofs = u[:, :neof]
    pcs = v[:, :neof]
    var = np.zeros([neof, 3])
    var[:, 0] = eigs[:neof]
    var[:, 1] = s[:neof] * 100
    var[:, 2] = var[:, 1] * np.sqrt(2/len(control))
    eofs = np.transpose(eofs)
    eofs = np.reshape(eofs, (neof, dim[1], dim[2]), order='F')
    return eofs, pcs, var


def eof_regress(pcs, eofn, control):
    '''
    Regresses original field on nth pcs time series

    Parameters
    ----------
    control: array
        original data to regress. time x nsigma x nlat
    eofn: int
        eof number to calculate
    pcs: array
        first principal component time series

    Returns
    -------
    eofn: array
        nth eof from field regressed onto pcsn
    '''
    from scipy import stats
    pcsn = pcs[:, eofn]
    pcsn = (pcsn - np.mean(pcsn)) / np.std(pcsn)
    eofn = np.zeros((np.ma.size(control, axis=1),
                    np.ma.size(control, axis=2)))
    for a in range(np.ma.size(control, axis=1)):
        for b in range(np.ma.size(control, axis=2)):
            eofn[a, b] = stats.linregress(pcsn, control[:, a, b])[0]
    return eofn  # this is not normalized!


def eof_response(response, eofn, lati, long, region=[-90, 90, 0, 360],
                 corr='no'):
    '''
    Projects response onto eofn

    Parameters
    ----------
    response: array
        data to project on to eofn
    eofn: array
        nth eof
    lat: array
        latitude of data
    lon: array
        longitude of data

    Returns
    -------
    projection: array
        response projected on eofn
    '''
    lat = lati[(region[0] <= lati) & (lati <= region[1])]
    lon = long[(region[2] <= long) & (long <= region[3])]
    eof_region = eofn.copy()[(region[0] <= lati) & (lati <= region[1]), :]
    eof_region = eof_region[:, (region[2] <= long) & (long <= region[3])]
    r_region = response.copy()[:, (region[0] <= lati) & (lati <= region[1]), :]
    r_region = r_region[:, :, (region[2] <= long) & (long <= region[3])]
    responsew = eofweight(r_region, lat, lon)
    eofn_w = eofweight(eof_region.copy(), lat, lon)
    projection = np.zeros(len(response))
    if corr == 'yes':
        for i in range(len(response)):
            projection[i] = (np.sum(responsew[i] * eofn_w) /
                             np.sqrt(np.sum(eofn_w*eofn_w) *
                                     np.sum(responsew[i]*responsew[i])))
    else:
        for i in range(len(response)):
            projection[i] = (np.sum(responsew[i] * eofn_w) /
                             (np.sum(eofn_w*eofn_w)))
    return projection


def eofweight(unweighted, lat, lon):
    '''
    Outputs weighted data for projections (i.e. sqrt cos lat)

    Parameters
    ----------
    unweighted: array
        unweighted data
    lat: array
        latitudes
    lon: array
        lon values

    Outputs
    -------
    weighted: array
        weighted data
    '''
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.sqrt(np.cos(meshlat.transpose() * np.pi/180))
    weighted = unweighted * meshlatweight
    return weighted
