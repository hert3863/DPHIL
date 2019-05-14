#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:45:23 2018

@author: bakerh
"""

import numpy as np
from netcdfread import ncread
##################
#### ALL ####
##################

def compare(items):
    '''
    Compares the control and perturbed ensembles
    and outputs a list of exp IDs that have completed
    for both ensembles, coupled with the patch number
    '''
    import glob
    import fnmatch
    # import lists of successful files for each exp
    controls = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/' + items +'/*')
    perturbs = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/' + items +'/*')
    # turn lists into just lists of exp IDs
    for i, item in enumerate(controls):
        controls[i] = controls[i][139+2*(len(items)-22):143+2*(len(items)-22)]
    for i, item in enumerate(perturbs):
        perturbs[i] = perturbs[i][139+2*(len(items)-22):143+2*(len(items)-22)]
    both = []
    # compare lists and add to dictionary if both exist
    for i, item in enumerate(controls):
        if fnmatch.filter(perturbs, controls[i]) != []:
            both.append(controls[i])
    return both


def dbase():
    from ANC import ANC
    lst = {}
    anc = ANC()
    anc.Start('b001')
    for i in range(1, 10001):
        lst[anc.Get()] = str(i)
        anc.Next()
    return lst


def calibration(recon, ncar):
    # computes calibration coefficients for reconstructions
    from scipy import stats
    ncar = np.reshape(ncar, (-1, np.shape(ncar)[-1]))
    recon = np.reshape(recon, (-1, np.shape(recon)[-1]))
    nind = np.shape(recon)[0]
    conv = np.zeros((nind, 2))
    for i in range(nind):
        conv[i, :] = stats.linregress(recon[i, 1:-5], ncar[i])[:2]
    return conv


def calibrate(recon, conv):
    # calibrates values
    oshape = np.shape(recon)
    recon = np.reshape(recon, (-1, np.shape(recon)[-1]))
    nind = np.shape(recon)[0]
    recon_corrected = np.zeros((nind, np.shape(recon)[1]))
    for i in range(nind):
        recon_corrected[i, :] = recon[i]*conv[i, 0] + conv[i, 1]
    recon_corrected = np.reshape(recon_corrected, oshape)
    return recon_corrected


##################
#### JET ####
##################


def lanczos(field, wi=51, co=0.0125, dim=1, hl='low'):
    def lweights(window, cutoff, hilo='low'):
        """
        Calculate weights for a low pass Lanczos filter
        and then applies to time series

        Parameters
        ----------

        window: int
            The length of the filter window.
        cutoff: float
            The cutoff frequency in inverse time steps.

        Returns
        -------
        w: array
            filter weights
        """
        order = ((window - 1) // 2) + 1
        nwts = 2 * order + 1
        w = np.zeros([nwts])
        n = nwts // 2
        w[n] = 2 * cutoff
        k = np.arange(1., n)
        sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
        firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
        w[n-1:0:-1] = firstfactor * sigma
        w[n+1:-1] = firstfactor * sigma
        w = w[1:-1]
        w = w / np.sum(w)
        if hilo == 'high':
            w = w*-1
            w[order-1] += 1
        return w

    l_filter = lweights(wi, co, hilo=hl)

    if dim == 1:
        filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1))
        for i in range(len(field)-len(l_filter)+1):
            filtered[i] = np.sum(field[i:i+len(l_filter)] * l_filter)
    if dim == 2:
        filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1,
                            np.ma.size(field, axis=1)))
        for i in range(len(field)-len(l_filter)+1):
            filtered[i] = np.sum((field[i:i+len(l_filter)].transpose() *
                                  l_filter).transpose(), axis=0)
    if dim == 3:
        filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1,
                            np.ma.size(field, axis=1),
                            np.ma.size(field, axis=2)))
        for i in range(len(field)-len(l_filter)+1):
            filtered[i] = np.sum((field[i:i+len(l_filter)].transpose((1, 2, 0)) *
                                  l_filter).transpose((2, 0, 1)), axis=0)
    return filtered


def jet_clim_indices():
    import glob
    from netcdfread import ncread
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/*')
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    ji = np.zeros((2, 2, len(a)))
    for i in range(len(a)):
        ucontrol = ncread(a[i], 'item15201_daily_mean')[:, 0, :]
        u_djf = ucontrol[300:450]
        u_jja = ucontrol[480:630]
        ji[0, 0, i], ji[0, 1, i] = jetind(u_djf, lat, lon)
        ji[1, 0, i], ji[1, 1, i] = jetind(u_jja, lat, lon)
        print(str(i))
    return ji


def jetind(u850, lat, lon):
    # need to pass in extra 30 days each side
    u850_96 = np.zeros((np.shape(u850)[0], 145, 192))
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    from scipy import interpolate
    for i in range(np.shape(u850)[0]):
        f = interpolate.interp2d(lon, lat[::-1], u850[i])
        u850_96[i] = f(lon_n96, lat_n96)
    u850_mean = np.mean(u850_96[:, 12:61, 160:], axis=2)
    u850_mean = lanczos(u850_mean, wi=61, co=1/10, dim=2)
    speeds = np.zeros((np.ma.size(u850_mean, axis=0)))
    lats = np.zeros((np.ma.size(u850_mean, axis=0)))
    lat_dist = lat_n96[12:61]
    for t in range(np.ma.size(u850_mean, axis=0)):
        speeds[t] = np.max(u850_mean[t, :])
        lats[t] = lat_dist[np.argmax(u850_mean[t, :])]
    speeds = np.mean(speeds)
    lats = np.mean(lats)
    return speeds, lats


def regress():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    both = compare('item15201_daily_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    ji = np.zeros((2, 2, len(both)))
    jic = np.zeros((2, 2, len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        u_djfc = ucontrol[300:450]
        u_jjac = ucontrol[480:630]
        u_djf = u[300:450]
        u_jja = u[480:630]
        jic[0, 0, i], jic[0, 1, i] = jetind(u_djfc, lat, lon)
        jic[1, 0, i], jic[1, 1, i] = jetind(u_jjac, lat, lon)
        ji[0, 0, i], ji[0, 1, i] = jetind(u_djf, lat, lon)
        ji[1, 0, i], ji[1, 1, i] = jetind(u_jja, lat, lon)
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    anom = anom[:, ::-1, :]
    jetind_d = ji - jic
    regmap_smoothed = np.zeros((2, 2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for a in range(2):
                for b in range(2):
                    regmap_smoothed[a, b, i+23, j] = np.cov(anom[:, i+23, j], jetind_d[a, b, :])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * np.cos(lat1[i+23]) * (np.pi/180)**2)
                    r = np.cov(anom[:, i+23, j], jetind_d[a, b, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], jetind_d[a, b, :])[1, 1]*np.cov(anom[:, i+23, j], jetind_d[a, b, :])[0, 0])
                    t = r*np.sqrt((len(both)-2)/(1-r**2))
                    p = 1 - stats.norm.cdf(np.abs(t))
                    sig = np.greater_equal(5, p*100*2).astype(int)
                    if sig == 1:
                        regmap_smoothed_sig[a, b, i+23, j] = regmap_smoothed[a, b, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def reconstruct_daily(gto, season='djf'):
    from scipy import interpolate
    from mpl_toolkits.basemap import shiftgrid
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 358.125, 192)
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    lon_sst_shift = np.linspace(.5, 359.5, 360)
    sindices = np.zeros((147*3))
    windices = np.zeros((147*3))
    if season == 'djf':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    if season == 'ndj':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [16+12*i, 17+12*i, 18+12*i]
            windices[3*i:3*(i+1)] = [10+12*i, 11+12*i, 12+12*i]
    if season == 'son':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [14+12*i, 15+12*i, 16+12*i]
            windices[3*i:3*(i+1)] = [8+12*i, 9+12*i, 10+12*i]
    if season == 'nd':
        months = 2
        sindices = np.zeros((147*2))
        windices = np.zeros((147*2))
        for i in range(147):
            sindices[2*i:2*(i+1)] = [16+12*i, 17+12*i]
            windices[2*i:2*(i+1)] = [10+12*i, 11+12*i]
    if season == 'dj':
        months = 2
        sindices = np.zeros((147*2))
        windices = np.zeros((147*2))
        for i in range(147):
            sindices[2*i:2*(i+1)] = [17+12*i, 18+12*i]
            windices[2*i:2*(i+1)] = [11+12*i, 12+12*i]
    if season == 'n':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 16+12*i
            windices[i] = 10+12*i
    if season == 'd':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 17+12*i
            windices[i] = 11+12*i
    if season == 'on':
        months = 2
        sindices = np.zeros((147*2))
        windices = np.zeros((147*2))
        for i in range(147):
            sindices[2*i:2*(i+1)] = [15+12*i, 16+12*i]
            windices[2*i:2*(i+1)] = [9+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    gto_interp = np.zeros((2, 2, 180, 360))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')
    sst_0 = np.where(sst<-999, 0, sst)
    for i in range(2):
        for j in range(2):
            f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto[i, j, :])
            gto_interp[i, j, :] = f(lon_sst_shift, lat_sst)
            gto_interp[i, j], lon1 = shiftgrid(181., gto_interp[i, j],
                                               lon_sst_shift, start=False)
    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2

    jet_recon_monthly = np.zeros((2, 2, months*147))
    jet_recon = np.zeros((2, 2, 147))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight
    for b in range(2):
        for y in range(147*months):
            jet_recon_monthly[0, b, y] = np.nansum(sst_w[y]*gto_interp[0, b])/20  # 2beta
            jet_recon_monthly[1, b, y] = np.nansum(sst_s[y]*gto_interp[1, b])/20  # 2beta

    if months == 1:
        for b in range(2):
            for y in range(147):
                jet_recon[0, b, y] = jet_recon_monthly[0, b, months*y]
                jet_recon[1, b, y] = jet_recon_monthly[1, b, months*y]
    else:
        for b in range(2):
            for y in range(147):
                jet_recon[0, b, y] = np.mean(jet_recon_monthly[0, b, months*y:months*(y+1)])
                jet_recon[1, b, y] = np.mean(jet_recon_monthly[1, b, months*y:months*(y+1)])

    return jet_recon

