#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:45:23 2018

@author: bakerh
"""

import numpy as np
from netcdfread import ncread


def main(dataset='ncep2'):
    def nao_region_ncep2(field):
        sindices = np.zeros((3*(np.shape(field)[0]//12-1)))
        windices = np.zeros((3*(np.shape(field)[0]//12-1)))
        for i in range((np.shape(field)[0]//12-1)):
            sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
        sindices = sindices.astype(int)
        windices = windices.astype(int)

        field1 = field[:, 4:29, 108:]
        field2 = field[:, 4:29, :17]
        field_r = np.concatenate((field1, field2), axis=2)
        field_w = field_r[windices]
        field_s = field_r[sindices]
        field_w = np.add.reduceat(field_w, range(0, len(field_w), 3))/3
        field_s = np.add.reduceat(field_s, range(0, len(field_s), 3))/3
        lat_r = np.arange(80, 17.5, -2.5)
        lon_r = np.arange(-90, 42.5, 2.5)
        return field_w, field_s, lat_r, lon_r

    def nao_region_ncar20c(field):
        sindices = np.zeros((3*(np.shape(field)[0]//12-1)))
        windices = np.zeros((3*(np.shape(field)[0]//12-1)))
        for i in range((np.shape(field)[0]//12-1)):
            sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
        sindices = sindices.astype(int)
        windices = windices.astype(int)

        field1 = field[:, 5:36, 135:]
        field2 = field[:, 5:36, :21]
        field_r = np.concatenate((field1, field2), axis=2)
        field_w = field_r[windices]
        field_s = field_r[sindices]
        field_w = np.add.reduceat(field_w, range(0, len(field_w), 3))/3
        field_s = np.add.reduceat(field_s, range(0, len(field_s), 3))/3
        lat_r = np.arange(80, 18, -2)
        lon_r = np.arange(-90, 42, 2)
        return field_w, field_s, lat_r, lon_r

    def nao_region_era20c(field):
        sindices = np.zeros((3*(np.shape(field)[0]//12-1)))
        windices = np.zeros((3*(np.shape(field)[0]//12-1)))
        for i in range((np.shape(field)[0]//12-1)):
            sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
        sindices = sindices.astype(int)
        windices = windices.astype(int)

        field1 = field[:, 10:71, 270:]
        field2 = field[:, 10:71, :41]
        field_r = np.concatenate((field1, field2), axis=2)
        field_w = field_r[windices]
        field_s = field_r[sindices]
        field_w = np.add.reduceat(field_w, range(0, len(field_w), 3))/3
        field_s = np.add.reduceat(field_s, range(0, len(field_s), 3))/3
        lat_r = np.arange(80, 19, -1)
        lon_r = np.arange(-90, 41, 1)
        return field_w, field_s, lat_r, lon_r

    from scipy import interpolate
    lat_n96 = np.arange(80, 18.75, -1.25)
    lon_n96 = np.linspace(-90, 39.375, 70)
    if dataset == 'ncep2':
        # NCEP2 MSLP 1979-2017 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'mslp')[:468] / 100
        mslp_w, mslp_s, lat_r, lon_r = nao_region_ncep2(mslp)
        mslp_w96 = np.zeros((np.shape(mslp_w)[0], 49, 70))
        mslp_s96 = np.zeros((np.shape(mslp_w)[0], 49, 70))
        for i in range(np.shape(mslp_w)[0]):
            f = interpolate.interp2d(lon_r, lat_r[::-1], mslp_w[i])
            mslp_w96[i] = f(lon_n96, lat_n96)
            f = interpolate.interp2d(lon_r, lat_r[::-1], mslp_s[i])
            mslp_s96[i] = f(lon_n96, lat_n96)
        eofs = np.zeros((2, 49, 70))
        var = np.zeros((2, 3))
        nao = np.zeros((2, len(mslp_w96)))
        eofs[0], var[0] = eof_clim(mslp_w96, lat_n96, lon_n96, 1)
        nao[0] = nao_indivyear(eofs[0], mslp_w96, lat_n96, lon_n96)
        eofs[1], var[1] = eof_clim(mslp_s96, lat_n96, lon_n96, 1)
        nao[1] = nao_indivyear(eofs[1], mslp_s96, lat_n96, lon_n96)

    if dataset == 'ncar20c':
        # NCAR-20thC MSLP 1871-2012 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc','prmsl') / 100
        mslp_w, mslp_s, lat_r, lon_r = nao_region_ncar20c(mslp)
        mslp_w96 = np.zeros((np.shape(mslp_w)[0], 49, 70))
        mslp_s96 = np.zeros((np.shape(mslp_w)[0], 49, 70))
        for i in range(np.shape(mslp_w)[0]):
            f = interpolate.interp2d(lon_r, lat_r[::-1], mslp_w[i])
            mslp_w96[i] = f(lon_n96, lat_n96)
            f = interpolate.interp2d(lon_r, lat_r[::-1], mslp_s[i])
            mslp_s96[i] = f(lon_n96, lat_n96)
        eofs = np.zeros((2, 49, 70))
        var = np.zeros((2, 3))
        nao = np.zeros((2, len(mslp_w96)))
        eofs[0], var[0] = eof_clim(mslp_w96, lat_n96, lon_n96, 1)
        nao[0] = nao_indivyear(eofs[0], mslp_w96, lat_n96, lon_n96)
        eofs[1], var[1] = eof_clim(mslp_s96, lat_n96, lon_n96, 1)
        nao[1] = nao_indivyear(eofs[1], mslp_s96, lat_n96, lon_n96)

    if dataset == 'era20c':
        # ERA20C MSLP 1900-2010 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc','msl') / 100
        mslp_w, mslp_s, lat_r, lon_r = nao_region_era20c(mslp)
        mslp_w96 = np.zeros((np.shape(mslp_w)[0], 49, 70))
        mslp_s96 = np.zeros((np.shape(mslp_w)[0], 49, 70))
        for i in range(np.shape(mslp_w)[0]):
            f = interpolate.interp2d(lon_r, lat_r[::-1], mslp_w[i])
            mslp_w96[i] = f(lon_n96, lat_n96)
            f = interpolate.interp2d(lon_r, lat_r[::-1], mslp_s[i])
            mslp_s96[i] = f(lon_n96, lat_n96)
        eofs = np.zeros((2, 49, 70))
        var = np.zeros((2, 3))
        nao = np.zeros((2, len(mslp_w96)))
        eofs[0], var[0] = eof_clim(mslp_w96, lat_n96, lon_n96, 1)
        nao[0] = nao_indivyear(eofs[0], mslp_w96, lat_n96, lon_n96)
        eofs[1], var[1] = eof_clim(mslp_s96, lat_n96, lon_n96, 1)
        nao[1] = nao_indivyear(eofs[1], mslp_s96, lat_n96, lon_n96)

    if dataset == 'hadam3p':
        import glob
        a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/*')
        mslp_w96 = np.zeros((len(a), 49, 70))
        mslp_s96 = np.zeros((len(a), 49, 70))
        for i, item in enumerate(a):
            mslp = ncread(item, 'item16222_monthly_mean')[:, 0, :] / 100
            mslp_w96[i], mslp_s96[i] = nao_region_hadam3p(mslp)
            print(str(i))
        eofs = np.zeros((2, 49, 70))
        var = np.zeros((2, 3))
        nao = np.zeros((2, len(mslp_w96)))
        eofs[0], var[0] = eof_clim(mslp_w96, lat_n96, lon_n96, 1)
        nao[0] = nao_indivyear(eofs[0], mslp_w96, lat_n96, lon_n96)
        eofs[1], var[1] = eof_clim(mslp_s96, lat_n96, lon_n96, 1)
        nao[1] = nao_indivyear(eofs[1], mslp_s96, lat_n96, lon_n96)
    return eofs, var, nao


def nao_region_hadam3p(field):
    sindices = [17, 18, 19]
    windices = [11, 12, 13]
    field1 = field[:, 8:57, 144:]
    field2 = field[:, 8:57, :22]
    field_r = np.concatenate((field1, field2), axis=2)
    field_w = np.mean(field_r[windices], axis=0)
    field_s = np.mean(field_r[sindices], axis=0)
    return field_w, field_s


def eof_corr(pat1, pat2, eof=0):
    lat = np.arange(80, 18.75, -1.25)
    lon = np.linspace(-90, 39.375, 70)
    pat1_w = eofweight(pat1, lat, lon)[eof]
    pat2_w = eofweight(pat2, lat, lon)[eof]
    crr = np.sum(pat1_w * pat2_w) / np.sqrt(np.sum(pat1_w * pat1_w) *
                                            np.sum(pat2_w * pat2_w))
    return crr


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


def eof_clim(clim, lat, lon, neof=2):
    def eof_svd(clim, lat, lon, neof):
        '''
        EOFs of data

        Parameters
        ----------
        clim: array
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
        field = clim.copy()
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
        var[:, 2] = var[:, 1] * np.sqrt(2/len(clim))
        eofs = np.transpose(eofs)
        eofs = np.reshape(eofs, (neof, dim[1], dim[2]), order='F')
        return eofs, pcs, var

    def eof_regress(pcs, eofn, clim):
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
        eofn = np.zeros((np.ma.size(clim, axis=1),
                        np.ma.size(clim, axis=2)))
        for a in range(np.ma.size(clim, axis=1)):
            for b in range(np.ma.size(clim, axis=2)):
                eofn[a, b] = stats.linregress(pcsn, clim[:, a, b])[0]
        return eofn  # this is not normalized!

    # compute eofs
    eofs_un, pcs, var = eof_svd(clim, lat, lon, neof)
    # regress control onto pcs to get sensible units
    eofs = np.zeros((neof, len(lat), len(lon)))
    for i in range(neof):
        eofs[i] = eof_regress(pcs, i, clim.copy())
    return eofs, var


def eof_response(response, eofn, lat, lon, corr='no'):
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
    eof_region = eofn.copy()
    r_region = response.copy()
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


def nao_indivyear(eof, indiv_years, lat, lon, ms='no'):
    # project responses onto eof
    proj = eof_response(indiv_years, eof, lat, lon)
    proj_m = np.mean(proj)
    proj_s = np.std(proj)
    proj = (proj-proj_m)/proj_s
    if ms == 'yes':
        return proj, proj_m, proj_s
    return proj


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


def regressnao(eof):
    from scipy import stats
    lat = np.arange(80, 18.75, -1.25)
    lon = np.linspace(-90, 39.375, 70)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    NAOinds = np.zeros((2, len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        mslp_control = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp_cw, mslp_cs = nao_region_hadam3p(mslp_control)
        mslp_w, mslp_s = nao_region_hadam3p(mslp)
        nao_cw = eof_response(np.expand_dims(mslp_cw, axis=0), eof[0], lat, lon)
        nao_cs = eof_response(np.expand_dims(mslp_cs, axis=0), eof[1], lat, lon)
        nao_w = eof_response(np.expand_dims(mslp_w, axis=0), eof[0], lat, lon)
        nao_s = eof_response(np.expand_dims(mslp_s, axis=0), eof[1], lat, lon)
        NAOinds[0, i] = nao_w - nao_cw
        NAOinds[1, i] = nao_s - nao_cs
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    anom = anom[:, ::-1, :]
    regmap_smoothed = np.zeros((2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for a in range(2):
                regmap_smoothed[a, i+23, j] = np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 1]*3/(4*6371**2 * 1.25 * 1.875 * (np.pi/180)**2)
                r = np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], NAOinds[a, :])[1, 1]*np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i+23, j] = regmap_smoothed[a, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def reconstruct_nao(gto, season='djf', rmask=np.ones((145, 192))):
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

    gto_interp = np.zeros((2, 180, 360))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')
    sst_0 = np.where(sst<-999, 0, sst)
    for i in range(2):
        f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto[i, :])
        gto_interp[i, :] = f(lon_sst_shift, lat_sst)
        gto_interp[i], lon1 = shiftgrid(181., gto_interp[i], lon_sst_shift,
                                        start=False)
    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2
    nao_recon_monthly = np.zeros((2, months*147))
    nao_recon = np.zeros((2, 147))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight
    '''
    OLD Rmask code
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    for y in range(36):
        nao_anom[0, y] = np.nansum(wdataSST[y]*gto[0]*(1-lsm1)*rmask)/20  # 2beta
        nao_anom[1, y] = np.nansum(sdataSST[y]*gto[1]*(1-lsm1)*rmask)/20
    '''
    for y in range(147*months):
        nao_recon_monthly[0, y] = np.nansum(sst_w[y]*gto_interp[0])/20  # 2beta
        nao_recon_monthly[1, y] = np.nansum(sst_s[y]*gto_interp[1])/20

    if months == 1:
        for y in range(147):
            nao_recon[0, y] = nao_recon_monthly[0, months*y]
            nao_recon[1, y] = nao_recon_monthly[1, months*y]
    else:
        for y in range(147):
            nao_recon[0, y] = np.mean(nao_recon_monthly[0, months*y:months*(y+1)])
            nao_recon[1, y] = np.mean(nao_recon_monthly[1, months*y:months*(y+1)])
    return nao_recon


def moving_corr(nao_ncep, nao_ncar, nao_era, nao_rec):
    ncep = np.zeros((2, 28))
    ncar = np.zeros((2, 131))
    era = np.zeros((2, 100))
    for s in range(2):
        nao_rec_s = (nao_rec[s] - np.mean(nao_rec[s])) / np.std(nao_rec[s])
        for i in range(28):
            ncep[s, i] = np.corrcoef(nao_ncep[s, i:i+11], nao_rec_s[109+i:109+i+11])[0, 1]
        for i in range(131):
            ncar[s, i] = np.corrcoef(nao_ncar[s, i:i+11], nao_rec_s[1+i:1+i+11])[0, 1]
        for i in range(100):
            era[s, i] = np.corrcoef(nao_era[s, i:i+11], nao_rec_s[30+i:30+i+11])[0, 1]
    return ncep, ncar, era
