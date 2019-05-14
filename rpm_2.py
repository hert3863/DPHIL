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


def regresspointwise():
    def win_sum(field):
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        field1 = field[:, :, 96:]
        field2 = field[:, :, :96]
        field_r = np.concatenate((field1, field2), axis=2)
        field_w = np.mean(field_r[windices], axis=0)
        field_s = np.mean(field_r[sindices], axis=0)
        return field_w, field_s

    from scipy import stats
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    pwinds = np.zeros((2, 145*192, len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        mslp_control = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp_cw, mslp_cs = win_sum(mslp_control)
        mslp_w, mslp_s = win_sum(mslp)
        pwinds[0, :, i] = mslp_w.flatten() - mslp_cw.flatten()
        pwinds[1, :, i] = mslp_s.flatten() - mslp_cs.flatten()
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    anom = anom[:, ::-1, :]
    regmap_smoothed = np.zeros((2, 145*192, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145*192, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for b in range(2):
                for a in range(145*192):
                    regmap_smoothed[b, a, i+23, j] = np.cov(anom[:, i+23, j], pwinds[b, a, :])[0, 1]*3/(4*6371**2 * 1.25 * 1.875 * (np.pi/180)**2)
                    r = np.cov(anom[:, i+23, j], pwinds[b, a, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], pwinds[b, a, :])[1, 1]*np.cov(anom[:, i+23, j], pwinds[b, a, :])[0, 0])
                    t = r*np.sqrt((len(both)-2)/(1-r**2))
                    p = 1 - stats.norm.cdf(np.abs(t))
                    sig = np.greater_equal(5, p*100*2).astype(int)
                    if sig == 1:
                        regmap_smoothed_sig[b, a, i+23, j] = regmap_smoothed[b, a, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


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


def mlr(nao, jet):
    import pandas as pd
    from statsmodels.formula.api import ols
    nao_st = np.zeros(np.shape(nao))
    jet_st = np.zeros(np.shape(jet))
    for i in range(2):
        nao_st[i] = (nao[i]-nao[i].mean())/np.std(nao[i])
        for j in range(2):
            jet_st[i, j] = (jet[i, j]-jet[i, j].mean())/np.std(jet[i, j])
    data1 = pd.DataFrame({'NAO': nao_st[0], 'Speed': jet_st[0, 0], 'Lat': jet_st[0, 1]})
    model1 = ols("NAO ~ Speed + Lat - 1", data1).fit()
    data2 = pd.DataFrame({'NAO': nao_st[1], 'Speed': jet_st[1, 0], 'Lat': jet_st[1, 1]})
    model2 = ols("NAO ~ Speed + Lat - 1", data2).fit()
    # print(model1.summary())
    # print(model._results.params)
    return model1._results.params, model2._results.params


def mlr_lanczos(nao, jet, nao_ncar, w=11):
    import pandas as pd
    from statsmodels.formula.api import ols
    w2 = int((w-1)*.5)
    nao_st = np.zeros(np.shape(nao))
    jet_st = np.zeros(np.shape(jet))
    for i in range(2):
        nao_st[i] = (nao[i]-nao[i].mean())/np.std(nao[i])
        for j in range(2):
            jet_st[i, j] = (jet[i, j]-jet[i, j].mean())/np.std(jet[i, j])
    jet_l = np.zeros((2, 2, np.shape(jet)[-1]-w+1))
    jet_h = np.zeros((2, 2, np.shape(jet)[-1]-w+1))
    for i in range(2):
        for j in range(2):
            jet_l[i, j] = lanczos(jet_st[i, j], w, 1/w)
            jet_h[i, j] = lanczos(jet_st[i, j], w, 1/w, hl='high')
    data1 = pd.DataFrame({'NAO': nao_st[0, w2:-w2], 'Speed_l': jet_l[0, 0], 'Lat_l': jet_l[0, 1], 'Speed_h': jet_h[0, 0], 'Lat_h': jet_h[0, 1]})
    model1 = ols("NAO ~ Speed_l + Speed_h + Lat_l +Lat_h - 1", data1).fit()
    data2 = pd.DataFrame({'NAO': nao_st[1, w2:-w2], 'Speed_l': jet_l[1, 0], 'Lat_l': jet_l[1, 1], 'Speed_h': jet_h[1, 0], 'Lat_h': jet_h[1, 1]})
    model2 = ols("NAO ~ Speed_l + Speed_h + Lat_l +Lat_h - 1", data2).fit()
    # print(model1.summary())
    # print(model._results.params)
    nao_w1111 = model1._results.params[0]*jet_l[0, 0] + model1._results.params[1]*jet_l[0, 1] + model1._results.params[2]*jet_h[0, 0] + model1._results.params[3]*jet_h[0, 1]
    nao_w0111 = model1._results.params[1]*jet_l[0, 1] + model1._results.params[2]*jet_h[0, 0] + model1._results.params[3]*jet_h[0, 1]
    nao_w1011 = model1._results.params[0]*jet_l[0, 0] + model1._results.params[2]*jet_h[0, 0] + model1._results.params[3]*jet_h[0, 1]
    nao_w1101 = model1._results.params[0]*jet_l[0, 0] + model1._results.params[1]*jet_l[0, 1] + model1._results.params[3]*jet_h[0, 1]
    nao_w1110 = model1._results.params[0]*jet_l[0, 0] + model1._results.params[1]*jet_l[0, 1] + model1._results.params[2]*jet_h[0, 0]
    nao_s1111 = model2._results.params[0]*jet_l[1, 0] + model2._results.params[1]*jet_l[1, 1] + model2._results.params[2]*jet_h[1, 0] + model2._results.params[3]*jet_h[1, 1]
    nao_s0111 = model2._results.params[1]*jet_l[1, 1] + model2._results.params[2]*jet_h[1, 0] + model2._results.params[3]*jet_h[1, 1]
    nao_s1011 = model2._results.params[0]*jet_l[1, 0] + model2._results.params[2]*jet_h[1, 0] + model2._results.params[3]*jet_h[1, 1]
    nao_s1101 = model2._results.params[0]*jet_l[1, 0] + model2._results.params[1]*jet_l[1, 1] + model2._results.params[3]*jet_h[1, 1]
    nao_s1110 = model2._results.params[0]*jet_l[1, 0] + model2._results.params[1]*jet_l[1, 1] + model2._results.params[2]*jet_h[1, 0]
    w_recons = np.zeros((2, 2))
    s_recons = np.zeros((2, 2))
    if w == 11:
        w3 = None
    else:
        w3 = -w2 + 5
    print(np.corrcoef(nao_w1111, nao_ncar[0, 4+w2-5:w3])[0, 1])
    print(np.corrcoef(nao_s1111, nao_ncar[1, 4+w2-5:w3])[0, 1])
    w_recons[0, 0] = np.corrcoef(nao_w0111, nao_ncar[0, 4+w2-5:w3])[0, 1]
    w_recons[0, 1] = np.corrcoef(nao_w1011, nao_ncar[0, 4+w2-5:w3])[0, 1]
    w_recons[1, 0] = np.corrcoef(nao_w1101, nao_ncar[0, 4+w2-5:w3])[0, 1]
    w_recons[1, 1] = np.corrcoef(nao_w1110, nao_ncar[0, 4+w2-5:w3])[0, 1]
    s_recons[0, 0] = np.corrcoef(nao_s0111, nao_ncar[1, 4+w2-5:w3])[0, 1]
    s_recons[0, 1] = np.corrcoef(nao_s1011, nao_ncar[1, 4+w2-5:w3])[0, 1]
    s_recons[1, 0] = np.corrcoef(nao_s1101, nao_ncar[1, 4+w2-5:w3])[0, 1]
    s_recons[1, 1] = np.corrcoef(nao_s1110, nao_ncar[1, 4+w2-5:w3])[0, 1]
    return w_recons, s_recons

##################
#### NAO ####
##################


def rean_eof(dataset='ncep2', n_eof=1):
    from scipy import interpolate
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    if dataset == 'ncep2':
        # NCEP2 MSLP 1979-2017 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'mslp')[:468] / 100
        lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'lon')
    elif dataset == 'ncar20c':
        # NCAR-20thC MSLP 1871-2012 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'prmsl') / 100
        lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'lon')
    elif dataset == 'era20c':
        # ERA20C MSLP 1900-2010 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc', 'msl') / 100
        lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc', 'latitude')
        lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc', 'longitude')

    mslp_96 = np.zeros((np.shape(mslp)[0], 145, 192))
    for i in range(np.shape(mslp)[0]):
        f = interpolate.interp2d(lon, lat[::-1], mslp[i])
        mslp_96[i] = f(lon_n96, lat_n96)

    sindices = np.zeros((3*(np.shape(mslp_96)[0]//12-1)))
    windices = np.zeros((3*(np.shape(mslp_96)[0]//12-1)))
    for i in range((np.shape(mslp_96)[0]//12-1)):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    field1 = mslp_96[:, :, 96:]
    field2 = mslp_96[:, :, :96]
    mslp_c = np.concatenate((field1, field2), axis=2)
    mslp_w = mslp_c[windices]
    mslp_s = mslp_c[sindices]
    mslp_w = np.add.reduceat(mslp_w, range(0, len(mslp_w), 3))/3
    mslp_s = np.add.reduceat(mslp_s, range(0, len(mslp_s), 3))/3

    eofs = np.zeros((2, 145, 192))
    var = np.zeros((2, 3))
    nao = np.zeros((2, len(mslp_w)))
    eofs[0], var[0] = eof_clim(mslp_w, lat_n96, lon_n96, n_eof)
    nao[0] = nao_indivyear(eofs[0, 8:57, 48:118], mslp_w[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    eofs[1], var[1] = eof_clim(mslp_s, lat_n96, lon_n96, n_eof)
    nao[1] = nao_indivyear(eofs[1, 8:57, 48:118], mslp_s[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    return eofs, var, nao


def had_eof(n_eof=1):
    import glob
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/*')
    mslp_w = np.zeros((len(a), 145, 192))
    mslp_s = np.zeros((len(a), 145, 192))
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    for i, item in enumerate(a):
        mslp_96 = ncread(item, 'item16222_monthly_mean')[:, 0, :] / 100
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        field1 = mslp_96[:, :, 96:]
        field2 = mslp_96[:, :, :96]
        mslp_c = np.concatenate((field1, field2), axis=2)
        mslp_w[i] = np.mean(mslp_c[windices], axis=0)
        mslp_s[i] = np.mean(mslp_c[sindices], axis=0)
        print(str(i))
    eofs = np.zeros((2, 145, 192))
    var = np.zeros((2, 3))
    nao = np.zeros((2, len(mslp_w)))
    eofs[0], var[0] = eof_clim(mslp_w, lat_n96, lon_n96, n_eof)
    nao[0] = nao_indivyear(eofs[0, 8:57, 48:118], mslp_w[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    eofs[1], var[1] = eof_clim(mslp_s, lat_n96, lon_n96, n_eof)
    nao[1] = nao_indivyear(eofs[1, 8:57, 48:118], mslp_s[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    return eofs, var, nao


def eof_corr(pat1, pat2):
    lat = np.arange(80, 18.75, -1.25)
    lon = np.linspace(-90, 39.375, 70)
    pat1_w = eofweight(pat1[8:57, 48:118], lat, lon)
    pat2_w = eofweight(pat2[8:57, 48:118], lat, lon)
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


def eof_clim(clim, lat, lon, neof=1, region=[8, 57, 48, 118]):
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
        neof -= 1
        from numpy import linalg
        field = clim.copy()[:, region[0]:region[1], region[2]:region[3]]
        dim = field.shape
        field = eofweight(field, lat[region[0]:region[1]],
                          lon[region[2]:region[3]])
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

        eofs = u[:, neof]
        pcs = v[:, neof]
        var = np.zeros(3)
        var[0] = eigs[neof]
        var[1] = s[neof] * 100
        var[2] = var[1] * np.sqrt(2/len(clim))
        eofs = np.transpose(eofs)
        eofs = np.reshape(eofs, (dim[1], dim[2]), order='F')
        return eofs, pcs, var

    def eof_regress(pcs, clim):
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
        pcsn = pcs[:]
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
    eofs = eof_regress(pcs, clim.copy())
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


def regressnao(eof):
    def win_sum(field):
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        field1 = field[:, :, 96:]
        field2 = field[:, :, :96]
        field_r = np.concatenate((field1, field2), axis=2)
        field_w = np.mean(field_r[windices], axis=0)
        field_s = np.mean(field_r[sindices], axis=0)
        return field_w, field_s

    from scipy import stats
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    NAOinds = np.zeros((2, len(both)))
    #NAOinds2 = np.zeros((2, len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        mslp_control = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp_cw, mslp_cs = win_sum(mslp_control)
        mslp_w, mslp_s = win_sum(mslp)
        nao_cw = eof_response(np.expand_dims(mslp_cw[8:57, 48:118], axis=0), eof[0, 8:57, 48:118], lat[8:57], lon[48:118])
        nao_cs = eof_response(np.expand_dims(mslp_cs[8:57, 48:118], axis=0), eof[1, 8:57, 48:118], lat[8:57], lon[48:118])
        nao_w = eof_response(np.expand_dims(mslp_w[8:57, 48:118], axis=0), eof[0, 8:57, 48:118], lat[8:57], lon[48:118])
        nao_s = eof_response(np.expand_dims(mslp_s[8:57, 48:118], axis=0), eof[1, 8:57, 48:118], lat[8:57], lon[48:118])
        NAOinds[0, i] = nao_w - nao_cw
        NAOinds[1, i] = nao_s - nao_cs

        #nao_cw2 = eof_response(np.expand_dims(mslp_cw[8:57, 48:118], axis=0), eof2[0, 8:57, 48:118], lat[8:57], lon[48:118])
        #nao_cs2 = eof_response(np.expand_dims(mslp_cs[8:57, 48:118], axis=0), eof2[1, 8:57, 48:118], lat[8:57], lon[48:118])
        #nao_w2 = eof_response(np.expand_dims(mslp_w[8:57, 48:118], axis=0), eof2[0, 8:57, 48:118], lat[8:57], lon[48:118])
        #nao_s2 = eof_response(np.expand_dims(mslp_s[8:57, 48:118], axis=0), eof2[1, 8:57, 48:118], lat[8:57], lon[48:118])
        #NAOinds2[0, i] = nao_w2 - nao_cw2
        #NAOinds2[1, i] = nao_s2 - nao_cs2
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    anom = anom[:, ::-1, :]
    regmap_smoothed = np.zeros((2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN

    #regmap_smoothed2 = np.zeros((2, 145, 192))
    #regmap_smoothed_sig2 = np.zeros((2, 145, 192))
    #regmap_smoothed_sig2[:] = np.NAN
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
                #regmap_smoothed2[a, i+23, j] = np.cov(anom[:, i+23, j], NAOinds2[a, :])[0, 1]*3/(4*6371**2 * 1.25 * 1.875 * (np.pi/180)**2)
                #r = np.cov(anom[:, i+23, j], NAOinds2[a, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], NAOinds2[a, :])[1, 1]*np.cov(anom[:, i+23, j], NAOinds2[a, :])[0, 0])
                #t = r*np.sqrt((len(both)-2)/(1-r**2))
                #p = 1 - stats.norm.cdf(np.abs(t))
                #sig = np.greater_equal(5, p*100*2).astype(int)
                #if sig == 1:
                #    regmap_smoothed_sig2[a, i+23, j] = regmap_smoothed2[a, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig#, regmap_smoothed2, regmap_smoothed_sig2


def reconstruct_nao(gto, season='djf', region='glob', detrend='no'):
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
    if season == 'ond':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [15+12*i, 16+12*i, 17+12*i]
            windices[3*i:3*(i+1)] = [9+12*i, 10+12*i, 11+12*i]
    if season == 'son':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [14+12*i, 15+12*i, 16+12*i]
            windices[3*i:3*(i+1)] = [8+12*i, 9+12*i, 10+12*i]
    if season == 'aso':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [13+12*i, 14+12*i, 15+12*i]
            windices[3*i:3*(i+1)] = [7+12*i, 8+12*i, 9+12*i]
    if season == 'jas':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [12+12*i, 13+12*i, 14+12*i]
            windices[3*i:3*(i+1)] = [6+12*i, 7+12*i, 8+12*i]
    if season == 'jja':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
            windices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    if season == 'mjj':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [10+12*i, 11+12*i, 12+12*i]
            windices[3*i:3*(i+1)] = [4+12*i, 5+12*i, 6+12*i]
    if season == 'amj':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [9+12*i, 10+12*i, 11+12*i]
            windices[3*i:3*(i+1)] = [3+12*i, 4+12*i, 5+12*i]
    if season == 'mam':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [8+12*i, 9+12*i, 10+12*i]
            windices[3*i:3*(i+1)] = [2+12*i, 3+12*i, 4+12*i]
    if season == 'fma':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [7+12*i, 8+12*i, 9+12*i]
            windices[3*i:3*(i+1)] = [1+12*i, 2+12*i, 3+12*i]
    if season == 'jfm':
        months = 3
        sindices = np.zeros((147*3))
        windices = np.zeros((147*3))
        for i in range(147):
            sindices[3*i:3*(i+1)] = [6+12*i, 7+12*i, 8+12*i]
            windices[3*i:3*(i+1)] = [0+12*i, 1+12*i, 2+12*i]
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
    if season == 'on':
        months = 2
        sindices = np.zeros((147*2))
        windices = np.zeros((147*2))
        for i in range(147):
            sindices[2*i:2*(i+1)] = [15+12*i, 16+12*i]
            windices[2*i:2*(i+1)] = [9+12*i, 11+12*i]
    if season == 'd':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 17+12*i
            windices[i] = 11+12*i
    if season == 'n':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 16+12*i
            windices[i] = 10+12*i
    if season == 'o':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 15+12*i
            windices[i] = 9+12*i
    if season == 's':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 14+12*i
            windices[i] = 8+12*i
    if season == 'a':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 13+12*i
            windices[i] = 7+12*i
    if season == 'j2':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 12+12*i
            windices[i] = 6+12*i
    if season == 'j':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 11+12*i
            windices[i] = 5+12*i
    if season == 'm':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 10+12*i
            windices[i] = 4+12*i
    if season == 'apr':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 9+12*i
            windices[i] = 3+12*i
    if season == 'mar':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 8+12*i
            windices[i] = 2+12*i
    if season == 'feb':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 7+12*i
            windices[i] = 1+12*i
    if season == 'jan':
        months = 1
        sindices = np.zeros((147))
        windices = np.zeros((147))
        for i in range(147):
            sindices[i] = 6+12*i
            windices[i] = 0+12*i
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    gto_interp = np.zeros((2, 180, 360))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')[:12*148]

    if detrend == 'yes':
        meshlat = np.zeros([np.ma.size(lon_sst), np.ma.size(lat_sst)])
        meshlat[:, :] = lat_sst
        meshlatweight = np.cos(meshlat.transpose() * np.pi/180)

        sst_glob_month = np.zeros((1776))
        for i in range(1776):
            test = np.where(sst[i] < -999, 0, 1)
            sst1 = np.where(sst[i] < -999, 0, sst[i])
            sst_glob_month[i] = np.sum((sst1*meshlatweight)[30:150])/(np.sum((meshlatweight*test)[30:150]))
        sst_glob_ts_sm = np.zeros((1776))
        for i in range(12):
            sst_month = sst_glob_month[i::12]
            sst_month1 = np.concatenate((sst_month[6:0:-1], sst_month, sst_month[-2:-8:-1]))
            sst_month_sm = lanczos(sst_month1, wi=13, co=1/16)
            sst_glob_ts_sm[i::12] = sst_month_sm
        sst = np.transpose(sst.transpose((1, 2, 0))-sst_glob_ts_sm, (2, 0, 1))

    sst_0 = np.where(sst < -999, 0, sst)
    gto_0 = np.where(np.isnan(gto), 0, gto)
    if region == 'atl':
        gto_0[:, :, 16:145] = 0 #90-0 for basin
    elif region == 'ind':
        gto_0[:, :, :16] = 0
        gto_0[:, :, 59:] = 0
    elif region == 'pac':
        gto_0[:, :, :59] = 0 #no equals fo basin
        gto_0[:, :, 145:] = 0 # 240E for basin
    for i in range(2):
        f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto_0[i, :])
        gto_interp[i, :] = f(lon_sst_shift, lat_sst)
        gto_interp[i], lon1 = shiftgrid(181., gto_interp[i], lon_sst_shift,
                                        start=False)

    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2
    nao_recon_monthly = np.zeros((2, months*147))
    nao_recon = np.zeros((2, 147))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight
    
    #gto_interp[:, 60:, :] = 0
    #gto_interp[1, :60, :] = 0
    #gto_interp[1, :, :280] = 0
    #gto_interp[1, :, 330:] = 0

    for y in range(147*months):
        nao_recon_monthly[0, y] = np.nansum(sst_w[y]*gto_interp[0])/20  # 2beta
        nao_recon_monthly[1, y] = np.nansum((sst_s[y])*gto_interp[1])/20

    if months == 1:
        for y in range(147):
            nao_recon[0, y] = nao_recon_monthly[0, months*y]
            nao_recon[1, y] = nao_recon_monthly[1, months*y]
    else:
        for y in range(147):
            nao_recon[0, y] = np.mean(nao_recon_monthly[0, months*y:months*(y+1)])
            nao_recon[1, y] = np.mean(nao_recon_monthly[1, months*y:months*(y+1)])
    #for i in range(2):
    #    nao_recon[i] = (nao_recon[i]-np.mean(nao_recon[i]))/np.std(nao_recon[i])
    '''
    sst_monthly = np.zeros((2, 147*3, 180, 360))
    for y in range(147*months):
        sst_monthly[0, y] = np.mean(sst_w[months*y:months*(y+1)])
        sst_monthly[1, y] = np.mean(sst_s[months*y:months*(y+1)])

    for y in range(147):
        nao_recon[0, y] = np.nansum(sst_monthly[0, y]*gto_interp[0])/20  # 2beta
        nao_recon[1, y] = np.nansum((sst_monthly[1, y])*gto_interp[1])/20
    '''
    return nao_recon


def laggedpredict(gto, nao_ncep, nao_ncar, nao_era):
    months = ['jan', 'feb', 'mar', 'apr', 'm', 'j', 'j2', 'a', 's', 'o', 'n',
              'd']
    regions = ['glob', 'atl', 'ind', 'pac']
    c_era = np.zeros((2, 4, 12))
    c_ncep = np.zeros((2, 4, 12))
    c_ncar = np.zeros((2, 4, 12))
    for i, month in enumerate(months):
        for j, region in enumerate(regions):
            nao_recon = reconstruct_nao(gto, season=month, region=region)
            for s in range(2):
                c_era[s, j, i] = np.corrcoef(nao_era[s], nao_recon[s, 30:140])[0, 1]
                c_ncep[s, j, i] = np.corrcoef(nao_ncep[s], nao_recon[s, 109:])[0, 1]
                c_ncar[s, j, i] = np.corrcoef(nao_ncar[s], nao_recon[s, 1:142])[0, 1]
        print(month)
    return c_ncep, c_ncar, c_era


def sstind():
    def moving_30(index):
        index30 = np.zeros((117))
        for i in range(117):
            index30[i] = np.mean(index[i:i+31])
        return index30

    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    sindices = np.zeros((147*3))
    for i in range(147):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
    sindices = sindices.astype(int)

    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180)

    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')
    sst_0 = np.where(sst < -999, 0, sst)
    mask = np.copy(sst)[0]
    mask = np.where(mask < -999, 0, 1)
    sst_1 = sst_0[sindices]
    sst_1 = np.mean(np.reshape(sst_1, (-1, 3, 180, 360)), axis=1)
    amo = (np.mean((sst_1*meshlatweight)[:, 20:90, 110:180], axis=(1,2)) / np.mean((meshlatweight*mask)[20:90, 110:180])-
           np.mean((sst_1*meshlatweight), axis=(1,2)) / np.mean((meshlatweight*mask)))
    n3 = np.mean((sst_1*meshlatweight)[:, 85:95, 30:90], axis=(1,2)) / np.mean((meshlatweight*mask)[85:95, 30:90])
    n34 = np.mean((sst_1*meshlatweight)[:, 85:95, 10:60], axis=(1,2)) / np.mean((meshlatweight*mask)[85:95, 10:60])
    amo_30 = moving_30(amo)
    n3_30 = moving_30(n3)
    n34_30 = moving_30(n34)
    amo_30 = (amo_30 - amo_30.mean())/np.std(amo_30)
    n3_30 = (n3_30 - n3_30.mean())/np.std(n3_30)
    n34_30 = (n34_30 - n34_30.mean())/np.std(n34_30)
    amo = (amo - amo.mean())/np.std(amo)
    n3 = (n3 - n3.mean())/np.std(n3)
    n34 = (n34 - n34.mean())/np.std(n34)
    return amo, n3, n34, amo_30, n3_30, n34_30


def montecarlotest(n):
    import scipy as sp
    from scipy import interpolate
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    # NCAR-20thC MSLP 1871-2012 inclusive
    mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'prmsl') / 100
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'lon')
    mslp_96 = np.zeros((np.shape(mslp)[0], 145, 192))
    for i in range(np.shape(mslp)[0]):
        f = interpolate.interp2d(lon, lat[::-1], mslp[i])
        mslp_96[i] = f(lon_n96, lat_n96)
    sindices = np.zeros((3*(np.shape(mslp_96)[0]//12)))
    for i in range((np.shape(mslp_96)[0]//12)):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    sindices = sindices.astype(int)
    field1 = mslp_96[:, :, 96:]
    field2 = mslp_96[:, :, :96]
    mslp_c = np.concatenate((field1, field2), axis=2)
    mslp_s = mslp_c[sindices]
    mslp_s = np.add.reduceat(mslp_s, range(0, len(mslp_s), 3))/3

    eof_control = eof_clim(mslp_s, lat_n96, lon_n96, 1)[0]
    if eof_control[32, 85] < 0:
        eof_control *= -1
    eof_control_w = eofweight(eof_control[8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    yr_ind = np.arange(0, 142)
    eof_42_12 = eof_clim(mslp_s[yr_ind[71:]], lat_n96, lon_n96, 1)[0]
    if eof_42_12[32, 85] < 0:
        eof_42_12 *= -1
    eof_42_12_w = eofweight(eof_42_12[8:57, 48:118], lat_n96[8:57], lon_n96[48:118])

    p = np.zeros(n+1)
    print('Starting loop')
    for j in range(n):
        yr_rand = np.random.permutation(yr_ind)
        yr_ind1 = yr_rand[:71]
        eof_n = eof_clim(mslp_s[yr_ind1], lat_n96, lon_n96, 1)[0]
        if eof_n[32, 85] < 0:
            eof_n *= -1
        eof_n_w = eofweight(eof_n[8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
        p[j] = (np.sum(eof_n_w * eof_control_w) /
                np.sqrt(np.sum(eof_n_w*eof_n_w) *
                        np.sum(eof_control_w*eof_control_w)))
        if np.remainder(j+1, 10) == 0:
            print(str(j+1) + '/' + str(n))
    p[-1] = (np.sum(eof_42_12_w * eof_control_w) /
             np.sqrt(np.sum(eof_42_12_w*eof_42_12_w) *
                     np.sum(eof_control_w*eof_control_w)))

    pos = sp.stats.percentileofscore(p[:-1], p[-1], kind='mean')
    sig = 1 - abs((pos-50)*2/100)
    # need sig less than 0.05 to be significant at 5% level
    return p, pos, sig


def coefficient_matrices():
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    psl = np.zeros((2, len(both), 145, 192))
    #ua = np.zeros((2, len(both), 144, 192))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        mslp_control = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        #u_control = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        #u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]

        windices = [11, 12, 13]
        sindices = [17, 18, 19]
        mslp_cs = np.mean(mslp_control[sindices], axis=0)
        mslp_s = np.mean(mslp[sindices], axis=0)
        mslp_cw = np.mean(mslp_control[windices], axis=0)
        mslp_w = np.mean(mslp[windices], axis=0)

        #u_cs = np.mean(u_control[510:600], axis=0)
        #u_s = np.mean(u[510:600], axis=0)
        #u_cw = np.mean(u_control[330:420], axis=0)
        #u_w = np.mean(u[330:420], axis=0)

        psl[0, i] = mslp_w - mslp_cw
        psl[1, i] = mslp_s - mslp_cs
        #ua[0, i] = u_w - u_cw
        #ua[1, i] = u_s - u_cs
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    anom = anom[:, ::-1, :]
    return anom, psl#, ua


def coefficient_regression(sst, z_sst, z_psl):
    # need to remove land points!!!
    y = np.reshape(sst, (27840))
    z = np.reshape(z_sst, (-1, 27840))
    x = np.zeros((5621, 27840))
    y = np.where(y < -999, 0, y)
    for i in range(5621):
        x[i] = np.where(y == 0, 0, z[i])
    bn = np.matmul(np.matmul(np.linalg.inv(np.matmul(x, x.transpose())), x), y)

    psl_recon = np.matmul(np.reshape(z_psl, (-1, 27840)).transpose(), bn)
    psl_recon = np.reshape(psl_recon, (145, 192))
    return psl_recon


def reconstruct_nao_gridskill(gto, nao_era):
    from scipy import interpolate
    from mpl_toolkits.basemap import shiftgrid
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 358.125, 192)
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    lon_sst_shift = np.linspace(.5, 359.5, 360)
    sindices = np.zeros((147*3))
    windices = np.zeros((147*3))

    months = 3
    sindices = np.zeros((147*3))
    windices = np.zeros((147*3))
    for i in range(147):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]

    sindices = sindices.astype(int)
    windices = windices.astype(int)

    gto_interp = np.zeros((2, 180, 360))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')
    sst_0 = np.where(sst < -999, 0, sst)
    gto_0 = np.where(np.isnan(gto), 0, gto)

    for i in range(2):
        f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto_0[i, :])
        gto_interp[i, :] = f(lon_sst_shift, lat_sst)
        gto_interp[i], lon1 = shiftgrid(181., gto_interp[i], lon_sst_shift,
                                        start=False)

    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2
    nao_recon_monthly = np.zeros((2, months*147, 180, 360))
    nao_recon = np.zeros((2, 147, 180, 360))
    nao_skill_map = np.zeros((2, 180, 360))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight

    for i in range(180):
        for j in range(360):
            for y in range(147*months):
                nao_recon_monthly[0, y, i, j] = np.nansum(sst_w[y, i, j]*gto_interp[0, i, j])/20  # 2beta
                nao_recon_monthly[1, y, i, j] = np.nansum((sst_s[y, i, j])*gto_interp[1, i, j])/20

            if months == 1:
                for y in range(147):
                    nao_recon[0, y, i, j] = nao_recon_monthly[0, months*y, i, j]
                    nao_recon[1, y, i, j] = nao_recon_monthly[1, months*y, i, j]
            else:
                for y in range(147):
                    nao_recon[0, y, i, j] = np.mean(nao_recon_monthly[0, months*y:months*(y+1), i, j])
                    nao_recon[1, y, i, j] = np.mean(nao_recon_monthly[1, months*y:months*(y+1), i, j])
            nao_skill_map[0, i, j] = np.corrcoef(nao_era[0, 24:94], nao_recon[0, 54:124, i, j])[0, 1]
            nao_skill_map[1, i, j] = np.corrcoef(nao_era[1, 24:94], nao_recon[0, 54:124, i, j])[0, 1]

    return nao_skill_map


def reconstruct_nao_regional(gto):
    from scipy import interpolate
    from mpl_toolkits.basemap import shiftgrid
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 358.125, 192)
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    lon_sst_shift = np.linspace(.5, 359.5, 360)
    sindices = np.zeros((147*3))
    windices = np.zeros((147*3))
    months = 3
    sindices = np.zeros((147*3))
    windices = np.zeros((147*3))
    for i in range(147):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]

    sindices = sindices.astype(int)
    windices = windices.astype(int)

    gto_interp = np.zeros((2, 180, 360))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')
    sst_0 = np.where(sst < -999, 0, sst)
    gto_0 = np.where(np.isnan(gto), 0, gto)

    for i in range(2):
        f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto_0[i, :])
        gto_interp[i, :] = f(lon_sst_shift, lat_sst)
        gto_interp[i], lon1 = shiftgrid(181., gto_interp[i], lon_sst_shift,
                                        start=False)

    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2
    nao_recon_monthly = np.zeros((8, months*147))
    nao_recon = np.zeros((8, 147))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight

    for y in range(147*months):
        nao_recon_monthly[0, y] = np.nansum(sst_s[y, 60:85, 225:280]*gto_interp[1, 60:85, 225:280])/20  # 2beta
        nao_recon_monthly[1, y] = np.nansum(sst_s[y, 85:110, 225:280]*gto_interp[1, 85:110, 225:280])/20
        nao_recon_monthly[2, y] = np.nansum(sst_s[y, 60:80, 280:330]*gto_interp[1, 60:80, 280:330])/20
        nao_recon_monthly[3, y] = np.nansum(sst_s[y, 80:110, 300:]*gto_interp[1, 80:110, 300:])/20 + np.nansum(sst_s[y, 80:110, :40]*gto_interp[1, 80:110, :40])/20
        nao_recon_monthly[4, y] = np.nansum(sst_s[y, 70:85, 60:90]*gto_interp[1, 70:85, 60:90])/20
        nao_recon_monthly[5, y] = np.nansum(sst_s[y, 80:110, 120:]*gto_interp[1, 80:110, 120:])/20
        nao_recon_monthly[7, y] = np.nansum(sst_s[y, 60:80, 330:]*gto_interp[1, 60:80, 330:])/20 + np.nansum(sst_s[y, 60:80, :40]*gto_interp[1, 60:80, :40])/20
        nao_recon_monthly[6, y] = np.nansum(sst_s[y, 40:60, 180:225]*gto_interp[1, 40:60, 180:225])/20

    for y in range(147):
        nao_recon[0, y] = np.mean(nao_recon_monthly[0, months*y:months*(y+1)])
        nao_recon[1, y] = np.mean(nao_recon_monthly[1, months*y:months*(y+1)])
        nao_recon[2, y] = np.mean(nao_recon_monthly[2, months*y:months*(y+1)])
        nao_recon[3, y] = np.mean(nao_recon_monthly[3, months*y:months*(y+1)])
        nao_recon[4, y] = np.mean(nao_recon_monthly[4, months*y:months*(y+1)])
        nao_recon[5, y] = np.mean(nao_recon_monthly[5, months*y:months*(y+1)])
        nao_recon[6, y] = np.mean(nao_recon_monthly[6, months*y:months*(y+1)])
        nao_recon[7, y] = np.mean(nao_recon_monthly[7, months*y:months*(y+1)])
    return nao_recon


##################
#### NAO STATION BASED 23rd October ####
##################

def rean_stationbased(dataset='ncep2'):
    from scipy import interpolate
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    if dataset == 'ncep2':
        # NCEP2 MSLP 1979-2017 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'mslp')[:468] / 100
        lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc', 'lon')
    elif dataset == 'ncar20c':
        # NCAR-20thC MSLP 1871-2012 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'prmsl')[:]/ 100
        lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/mslp/prmsl.mon.mean.nc', 'lon')
    elif dataset == 'era20c':
        # ERA20C MSLP 1900-2010 inclusive
        mslp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc', 'msl') / 100
        lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc', 'latitude')
        lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/ERA20C_mslp_daily_1900_2010.nc', 'longitude')

    mslp_96 = np.zeros((np.shape(mslp)[0], 145, 192))
    for i in range(np.shape(mslp)[0]):
        f = interpolate.interp2d(lon, lat[::-1], mslp[i])
        mslp_96[i] = f(lon_n96, lat_n96)

    sindices = np.zeros((3*(np.shape(mslp_96)[0]//12-1)))
    windices = np.zeros((3*(np.shape(mslp_96)[0]//12-1)))
    for i in range((np.shape(mslp_96)[0]//12-1)):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    mslp_s = mslp_96[sindices]
    mslp_w = mslp_96[windices]
    mslp_s = np.add.reduceat(mslp_s, range(0, len(mslp_s), 3))/3
    mslp_w = np.add.reduceat(mslp_w, range(0, len(mslp_w), 3))/3

    nao_s = mslp_s[:, 28, 191] - mslp_s[:, 18, 168]
    nao_s = (nao_s - nao_s.mean())/np.std(nao_s)
    nao_w = mslp_w[:, 41, 187] - mslp_w[:, 20, 180]
    nao_w = (nao_w - nao_w.mean())/np.std(nao_w)
    nao = np.stack((nao_w, nao_s))
    return nao


def regressnao_stationbased():
    from scipy import stats
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    NAOinds = np.zeros((2, len(both)))
    nao_c = np.zeros((2, len(both)))
    nao = np.zeros((2, len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        mslp_control = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        mslp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        windices = [11, 12, 13]
        sindices = [17, 18, 19]
        mslp_cs = np.mean(mslp_control[sindices], axis=0)
        mslp_s = np.mean(mslp[sindices], axis=0)
        mslp_cw = np.mean(mslp_control[windices], axis=0)
        mslp_w = np.mean(mslp[windices], axis=0)
        nao_c[1, i] = mslp_cs[36, 176] - mslp_cs[18, 168]  #mslp_cs[28, 191] - mslp_cs[18, 168]
        nao_c[0, i] = mslp_cw[41, 187] - mslp_cw[20, 180]
        nao[1, i] = mslp_s[36, 176] - mslp_s[18, 168]  #mslp_s[28, 191] - mslp_s[18, 168]
        nao[0, i] = mslp_w[41, 187] - mslp_w[20, 180]

        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    nao[0] = (nao[0] - nao_c[0].mean())/np.std(nao_c[0])
    nao[1] = (nao[1] - nao_c[1].mean())/np.std(nao_c[1])
    nao_c[0] = (nao_c[0] - nao_c[0].mean())/np.std(nao_c[0])
    nao_c[1] = (nao_c[1] - nao_c[1].mean())/np.std(nao_c[1])
    NAOinds = nao - nao_c
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


def jetpos():
    import glob
    from netcdfread import ncread
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/*')
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    ji = np.zeros((2, 192))
    for j in range(len(a)):
        ucontrol = ncread(a[j], 'item15201_daily_mean')[:, 0, :]
        u_djf = ucontrol[330:420]
        u_jja = ucontrol[510:600]
        lat_n96 = np.linspace(90, -90, 145)
        lon_n96 = np.arange(0, 360, 360/192)
        u_djf_96 = np.zeros((np.shape(u_djf)[0], 145, 192))
        u_jja_96 = np.zeros((np.shape(u_djf)[0], 145, 192))
        lat_dist = lat_n96[12:61]
        from scipy import interpolate
        for i in range(np.shape(u_djf)[0]):
            f = interpolate.interp2d(lon, lat[::-1], u_djf[i])
            u_djf_96[i] = f(lon_n96, lat_n96)
            f = interpolate.interp2d(lon, lat[::-1], u_jja[i])
            u_jja_96[i] = f(lon_n96, lat_n96)
        u_jja_96 = np.mean(u_jja_96[:, 12:61], axis=0)
        u_djf_96 = np.mean(u_djf_96[:, 12:61], axis=0)
        ji[1, :] += lat_dist[np.argmax(u_jja_96, axis=0)]
        ji[0, :] += lat_dist[np.argmax(u_djf_96, axis=0)]
        print(str(j))
    ji /= len(a)
    return ji


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


def stpos():
    from netcdfread import ncread

    def comparesame():
        '''
        Compares the control and perturbed ensembles
        and outputs a list of exp IDs that have completed
        for both ensembles, coupled with the patch number
        '''
        import glob
        import fnmatch
        # import lists of successful files for each exp
        controls = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_daily_mean/*')
        perturbs = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16203_daily_mean/*')
        # turn lists into just lists of exp IDs
        for i, item in enumerate(controls):
            controls[i] = controls[i][139+2*(20-22):143+2*(20-22)]
        for i, item in enumerate(perturbs):
            perturbs[i] = perturbs[i][139+2*(20-22):143+2*(20-22)]
        both = []
        # compare lists and add to dictionary if both exist
        for i, item in enumerate(controls):
            if fnmatch.filter(perturbs, controls[i]) != []:
                both.append(controls[i])
        return both

    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    both = comparesame()
    sti = np.zeros((2, 192))
    for j in range(len(both)):
        v = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_daily_mean/item15202_daily_mean_' + both[j] + '_1999-01_2000-12.nc', 'item15202_daily_mean')[:, 0, :]
        t = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16203_daily_mean/item16203_daily_mean_' + both[j] + '_1999-01_2000-12.nc', 'item16203_daily_mean')[:, 0, :]
        v_djf = v[300:450]
        v_jja = v[480:630]
        t_djf_96 = t[300:450]
        t_jja_96 = t[480:630]
        lat_n96 = np.linspace(90, -90, 145)
        lon_n96 = np.arange(0, 360, 360/192)
        v_djf_96 = np.zeros((np.shape(v_djf)[0], 145, 192))
        v_jja_96 = np.zeros((np.shape(v_djf)[0], 145, 192))
        lat_dist = lat_n96[12:61]
        from scipy import interpolate
        for i in range(np.shape(v_djf)[0]):
            f = interpolate.interp2d(lon, lat[::-1], v_djf[i])
            v_djf_96[i] = f(lon_n96, lat_n96)
            f = interpolate.interp2d(lon, lat[::-1], v_jja[i])
            v_jja_96[i] = f(lon_n96, lat_n96)
        v_jja_96 = lanczos(v_jja_96, 61, 1/6, dim=3, hl='high')
        v_djf_96 = lanczos(v_djf_96, 61, 1/6, dim=3, hl='high')
        t_jja_96 = lanczos(t_jja_96, 61, 1/6, dim=3, hl='high')
        t_djf_96 = lanczos(t_djf_96, 61, 1/6, dim=3, hl='high')
        vt_jja = (v_jja_96 * t_jja_96).mean(axis=0)
        vt_djf = (v_djf_96 * t_djf_96).mean(axis=0)
        sti[1, :] += lat_dist[np.argmax(vt_jja[12:61], axis=0)]
        sti[0, :] += lat_dist[np.argmax(vt_djf[12:61], axis=0)]
        print(str(j))
    sti /= len(both)
    return sti


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


def jetind_rean():
    import glob
    from netcdfread import ncread
    jet_ncep = np.zeros((2, 2, 38))
    jet_ncar = np.zeros((2, 2, 141))
    jet_era = np.zeros((2, 2, 110))

    a = glob.glob('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/*')
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/uwnd.1979.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/uwnd.1979.nc','lon')
    for i in range(len(a)-1):
        u_jf = ncread(a[i+1], 'uwnd')[:89, 2, :]
        u_jja = ncread(a[i+1], 'uwnd')[-244:-92, 2, :]
        u_d = ncread(a[i], 'uwnd')[-61:, 2, :]
        u_djf = np.concatenate((u_d, u_jf), axis=0)
        jet_ncep[0, 0, i], jet_ncep[0, 1, i] = jetind(u_djf, lat, lon)
        jet_ncep[1, 0, i], jet_ncep[1, 1, i] = jetind(u_jja, lat, lon)
        print('ncep' + str(i))

    a = glob.glob('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/u850/*')
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/u850/uwnd.1871.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/u850/uwnd.1871.nc','lon')
    for i in range(len(a)-1):
        u_jf = ncread(a[i+1], 'uwnd')[:89, 3, :]
        u_jja = ncread(a[i+1], 'uwnd')[-244:-92, 3, :]
        u_d = ncread(a[i], 'uwnd')[-61:, 3, :]
        u_djf = np.concatenate((u_d, u_jf), axis=0)
        jet_ncar[0, 0, i], jet_ncar[0, 1, i] = jetind(u_djf, lat, lon)
        jet_ncar[1, 0, i], jet_ncar[1, 1, i] = jetind(u_jja, lat, lon)
        print('ncar' + str(i))

    a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/u850/*'))
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/u850/ERA20C_U850_dailymeaned_1900.nc','g0_lat_1')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/ERA20C/u850/ERA20C_U850_dailymeaned_1900.nc','g0_lon_2')
    for i in range(len(a)-1):
        u_jf = ncread(a[i+1], 'U_GDS0_ISBL')[:89, :]
        u_jja = ncread(a[i+1], 'U_GDS0_ISBL')[-244:-92, :]
        u_d = ncread(a[i], 'U_GDS0_ISBL')[-61:, :]
        u_djf = np.concatenate((u_d, u_jf), axis=0)
        jet_era[0, 0, i], jet_era[0, 1, i] = jetind(u_djf, lat, lon)
        jet_era[1, 0, i], jet_era[1, 1, i] = jetind(u_jja, lat, lon)
        print('era' + str(i))

    return jet_ncep, jet_ncar, jet_era


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

##################
#### SCHALLER ####
##################

def regress_je():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    both = compare('item15201_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    jec = np.zeros((len(both)))
    je = np.zeros((len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        u_djfc = ucontrol[11:14].mean(axis=0)
        u_djf = u[11:14].mean(axis=0)
        jec[i] = 0.5*(np.mean((u_djfc*meshlatweight)[24:36, 184:]) / np.mean(meshlatweight[24:36, 184:]) +
                      np.mean((u_djfc*meshlatweight)[24:36, :8]) / np.mean(meshlatweight[24:36, :8]))
        je[i] = 0.5*(np.mean((u_djf*meshlatweight)[24:36, 184:]) / np.mean(meshlatweight[24:36, 184:]) +
                      np.mean((u_djf*meshlatweight)[24:36, :8]) / np.mean(meshlatweight[24:36, :8]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    anom = anom[:, ::-1, :]
    je_d = je - jec
    regmap_smoothed = np.zeros((145, 192))
    regmap_smoothed_sig = np.zeros((145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            regmap_smoothed[i+23, j] = np.cov(anom[:, i+23, j], je_d[:])[0, 1]*3/(4 * 6371**2 * 1.25 * 1.875 * np.cos(lat1[i+23]) * (np.pi/180)**2)
            r = np.cov(anom[:, i+23, j], je_d[:])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], je_d[:])[1, 1]*np.cov(anom[:, i+23, j], je_d[:])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                regmap_smoothed_sig[i+23, j] = regmap_smoothed[i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def sst_schaller(gtojet):
    import glob
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/schaller_ssts/*')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    meshlatweight_gto = np.cos(meshlat.transpose() * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2
    gradient = np.zeros(11)
    js = np.zeros(11)
    for i, pat in enumerate(a):
        sst = ncread(pat, 'tos')[0, ::-1]
        mask = np.copy(sst)
        mask = np.where(mask > 1000, 0, 1)
        sst_n = np.mean((sst*meshlatweight)[16:33, 171:]) / np.mean((meshlatweight*mask)[16:33, 171:])
        sst_s = np.mean((sst*meshlatweight)[32:49, 171:]) / np.mean((meshlatweight*mask)[32:49, 171:])
        gradient[i] = sst_s - sst_n
        js[i] = np.nansum((sst*gtojet*meshlatweight_gto)[24:61, 160:])/20
    return gradient, js


def schaller_na(gtojet):
    import glob
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/schaller_ssts/*')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight_gto = np.cos(meshlat.transpose() * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2
    js_na = np.zeros(11)
    js = np.zeros(11)
    for i, pat in enumerate(a):
        sst = ncread(pat, 'tos')[0, ::-1]
        js_na[i] = np.nansum((sst*gtojet*meshlatweight_gto)[24:61, 160:])/20
        js[i] = np.nansum((sst*gtojet*meshlatweight_gto))/20
    return js_na, js


def schaller_na2(gtojet, jete_recon_corr):
    # same as above but looking at anomalies from background Jan 2014 SST state
    import glob
    from netcdfread import ncread
    import scipy.stats as stats
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    sst_hist = ncread('/network/aopp/hera/mad/bakerh/RPM/happissts/ALLclim_ancil_146months_OSTIA_sst_2004-12-01_2017-01-30.nc', 'temp')[:864, 0]
    sst_hist = np.mean(sst_hist.reshape((-1, 6, 145, 192)), axis=1)
    sst_2014 = sst_hist[109]
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/schaller_ssts/*')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight_gto = np.cos(meshlat.transpose() * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2
    js_na = np.zeros(12)
    js_in = np.zeros(12)
    js_pa = np.zeros(12)
    js = np.zeros(12)
    sst_all = {}
    # ensure all sst patterns have the same mask
    for i, pat in enumerate(a):
        sst_all[pat] = ncread(pat, 'tos')[0, ::-1]
    for i, pat in enumerate(a):
        sst = sst_all[pat] * 0
        for pa in a:
            sst += sst_all[pa] * 0
        js_na[i] = np.nansum(((sst_2014-sst_all[pat]+sst-273.15)*gtojet*meshlatweight_gto)[24:61, 160:])/20
        js_in[i] = np.nansum(((sst_2014-sst_all[pat]+sst-273.15)*gtojet*meshlatweight_gto)[60:97, 16:57])/20
        js_pa[i] = np.nansum(((sst_2014-sst_all[pat]+sst-273.15)*gtojet*meshlatweight_gto)[48:89, 57:129])/20
        js[i] = np.nansum(((sst_2014-sst_all[pat]+sst-273.15)*gtojet*meshlatweight_gto))/20
    js_na[11] = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto)[24:61, 160:])/20
    js_in[11] = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto)[60:97, 16:57])/20
    js_pa[11] = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto)[48:89, 57:129])/20
    js[11] = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto))/20

    sst_djf = np.mean(np.reshape(sst_hist, (-1, 3, 145, 192)), axis=1)[0::4]
    je_happi = np.zeros(12)
    for i in range(12):
        je_happi[i] = np.nansum(((sst_djf[i]+sst-273.15)*gtojet*meshlatweight_gto))/20
    conv = stats.linregress(je_happi, jete_recon_corr[-13: -1]-jete_recon_corr.mean())
    js_na1 = js_na / js
    js_in1 = js_in / js
    js_pa1 = js_pa / js
    js1 = js * conv[0] + conv[1]
    js_na1 *= js1
    js_in1 *= js1
    js_pa1 *= js1
    return js_na1, js_in1, js_pa1, js1


def schaller_na3(gtojet, jete_recon_corr):
    # same as above but looking at anomalies from 0 then
    # subtracting from total of jan 2014
    import glob
    from netcdfread import ncread
    import scipy.stats as stats
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    sst_hist = ncread('/network/aopp/hera/mad/bakerh/RPM/happissts/ALLclim_ancil_146months_OSTIA_sst_2004-12-01_2017-01-30.nc', 'temp')[:864, 0]
    sst_hist = np.mean(sst_hist.reshape((-1, 6, 145, 192)), axis=1)
    sst_2014 = sst_hist[109]
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/schaller_ssts/*')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight_gto = np.cos(meshlat.transpose() * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2
    js_na = np.zeros(11)
    js_in = np.zeros(11)
    js_pa = np.zeros(11)
    js = np.zeros(11)
    sst_all = {}
    # ensure all sst patterns have the same mask
    for i, pat in enumerate(a):
        sst_all[pat] = ncread(pat, 'tos')[0, ::-1]
    for i, pat in enumerate(a):
        sst = sst_all[pat] * 0
        for pa in a:
            sst += sst_all[pa] * 0
        js_na[i] = np.nansum(((sst_2014*0+sst_all[pat]+sst)*gtojet*meshlatweight_gto)[24:61, 160:])/20
        js_in[i] = np.nansum(((sst_2014*0+sst_all[pat]+sst)*gtojet*meshlatweight_gto)[60:97, 16:57])/20
        js_pa[i] = np.nansum(((sst_2014*0+sst_all[pat]+sst)*gtojet*meshlatweight_gto)[48:89, 57:129])/20
        js[i] = np.nansum(((sst_2014-sst_all[pat]+sst-273.15)*gtojet*meshlatweight_gto))/20
    js_na_hist = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto)[24:61, 160:])/20
    js_in_hist = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto)[60:97, 16:57])/20
    js_pa_hist = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto)[48:89, 57:129])/20
    js_hist = np.nansum(((sst_2014+sst-273.15)*gtojet*meshlatweight_gto))/20

    js_na = js_hist - js_na
    js_in = js_hist - js_in
    js_pa = js_hist - js_pa

    sst_djf = np.mean(np.reshape(sst_hist, (-1, 3, 145, 192)), axis=1)[0::4]
    je_happi = np.zeros(12)
    for i in range(12):
        je_happi[i] = np.nansum(((sst_djf[i]+sst-273.15)*gtojet*meshlatweight_gto))/20
    conv = stats.linregress(je_happi, jete_recon_corr[-13: -1]-jete_recon_corr.mean())

    js_na_hist = js_na_hist / js_hist
    js_in_hist = js_in_hist / js_hist
    js_pa_hist = js_pa_hist / js_hist

    js1 = js * conv[0] + conv[1]
    js_hist = js_hist * conv[0] + conv[1]
    js_na1 = js_na * conv[0] + conv[1]
    js_in1 = js_in * conv[0] + conv[1]
    js_pa1 = js_pa * conv[0] + conv[1]

    js_na_hist *= js_hist
    js_in_hist *= js_hist
    js_pa_hist *= js_hist

    js_ind = np.hstack((js_in1, js_hist-js_in_hist))
    js_atl = np.hstack((js_na1, js_hist-js_na_hist))
    js_pac = np.hstack((js_pa1, js_hist-js_pa_hist))
    js_all = np.hstack((js1, js_hist))


    return js_all, js_atl, js_ind, js_pac


def jetextind(u200, lat, lon):
    # pass in u200
    u_mean = u200.mean(axis=0)
    u_96 = np.zeros((145, 192))
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    meshlat = np.zeros([np.ma.size(lon_n96), np.ma.size(lat_n96)])
    meshlat[:, :] = lat_n96
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    from scipy import interpolate
    f = interpolate.interp2d(lon, lat[::-1], u_mean)
    u_96 = f(lon_n96, lat_n96)
    je = 0.5*(np.mean((u_96*meshlatweight)[24:36, 184:]) / np.mean(meshlatweight[24:36, 184:]) +
              np.mean((u_96*meshlatweight)[24:36, :8]) / np.mean(meshlatweight[24:36, :8]))
    return je


def jetextind_rean():
    import glob
    from netcdfread import ncread
    jete_ncar = np.zeros((141))
    jete_ncep = np.zeros((38))

    a = glob.glob('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/u850/*')
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/u850/uwnd.1871.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCAR-20thC/u850/uwnd.1871.nc','lon')
    for i in range(len(a)-1):
        u_jf = ncread(a[i+1], 'uwnd')[:59, 16, :]
        u_d = ncread(a[i], 'uwnd')[-31:, 16, :]
        u_djf = np.concatenate((u_d, u_jf), axis=0)
        jete_ncar[i] = jetextind(u_djf, lat, lon)
        print('ncar' + str(i))

    a = glob.glob('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/*')
    lat = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/uwnd.1979.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/uwnd.1979.nc','lon')
    for i in range(len(a)-1):
        u_jf = ncread(a[i+1], 'uwnd')[:59, 9, :]
        u_d = ncread(a[i], 'uwnd')[-31:, 9, :]
        u_djf = np.concatenate((u_d, u_jf), axis=0)
        jete_ncep[i] = jetextind(u_djf, lat, lon)
        print('ncep' + str(i))
    return jete_ncep, jete_ncar


def reconstruct_dailyje(gto, season='djf'):
    from scipy import interpolate
    from mpl_toolkits.basemap import shiftgrid
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 358.125, 192)
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    lon_sst_shift = np.linspace(.5, 359.5, 360)
    windices = np.zeros((147*3))
    if season == 'djf':
        months = 3
        windices = np.zeros((147*3))
        for i in range(147):
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]

    windices = windices.astype(int)

    gto_interp = np.zeros((180, 360))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')
    sst_0 = np.where(sst <- 999, 0, sst)
    f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto[:])
    gto_interp[:] = f(lon_sst_shift, lat_sst)
    gto_interp, lon1 = shiftgrid(181., gto_interp, lon_sst_shift, start=False)
    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2

    jet_recon_monthly = np.zeros((months*147))
    jet_recon = np.zeros((147))

    sst_w = sst_0[windices]
    sst_w = sst_w * meshlatweight
    for y in range(147*months):
        jet_recon_monthly[y] = np.nansum((sst_w[y]*gto_interp))/20  # 2beta

    if months == 1:
        for y in range(147):
            jet_recon[y] = jet_recon_monthly[months*y]
    else:
        for y in range(147):
            jet_recon[y] = np.mean(jet_recon_monthly[months*y:months*(y+1)])

    return jet_recon

##################
#### FUTURE ####
##################

def happi_recon(gto, season='djf', region='glob'):
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 358.125, 192)
    sindices = np.zeros((12*3))
    windices = np.zeros((12*3))
    if season == 'djf':
        months = 3
        sindices = np.zeros((12*3))
        windices = np.zeros((12*3))
        for i in range(12):
            sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
            windices[3*i:3*(i+1)] = [0+12*i, 1+12*i, 2+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    sst_hist = ncread('/network/aopp/hera/mad/bakerh/RPM/happissts/ALLclim_ancil_146months_OSTIA_sst_2004-12-01_2017-01-30.nc', 'temp')[:870, 0]
    sst_hist = np.mean(sst_hist.reshape((-1, 6, 145, 192)), axis=1)
    sst_15 = ncread('/network/aopp/hera/mad/bakerh/RPM/happissts/HAPPI15C_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30.nc', 'temp')[:870, 0]
    sst_15 = np.mean(sst_15.reshape((-1, 6, 145, 192)), axis=1)
    sst_20 = ncread('/network/aopp/hera/mad/bakerh/RPM/happissts/HAPPI20C_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30.nc', 'temp')[:870, 0]
    sst_20 = np.mean(sst_20.reshape((-1, 6, 145, 192)), axis=1)
    sst_30 = ncread('/network/aopp/hera/mad/bakerh/RPM/happissts/HAPPI30_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30.nc', 'temp')[:870, 0]
    sst_30 = np.mean(sst_30.reshape((-1, 6, 145, 192)), axis=1)

    sst_hist0 = np.where(sst_hist < -999, 0, sst_hist)
    sst_150 = np.where(sst_15 < -999, 0, sst_15)
    sst_200 = np.where(sst_20 < -999, 0, sst_20)
    sst_300 = np.where(sst_30 < -999, 0, sst_30)

    gto_0 = np.where(np.isnan(gto), 0, gto)
    if region == 'atl':
        gto_0[:, :, lon_n96 < 270] = 0
    elif region == 'ind':
        gto_0[:, :, lon_n96 < 30] = 0
        gto_0[:, :, lon_n96 > 110] = 0
    elif region == 'pac':
        gto_0[:, :, lon_n96 < 110] = 0
        gto_0[:, :, lon_n96 > 240] = 0

    meshlat = np.zeros([145, 192])
    meshlat[:, :] = np.expand_dims(lat_n96, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2

    years = [sst_hist0, sst_150, sst_200, sst_300]
    nao_recons = np.zeros((4, 2, 12))
    for s, sst_0 in enumerate(years):
        nao_recon_monthly = np.zeros((2, months*12))
        nao_recon = np.zeros((2, 12))

        sst_w = sst_0[windices]
        sst_s = sst_0[sindices]
        sst_w = sst_w * meshlatweight
        sst_s = sst_s * meshlatweight

        for y in range(12*months):
            nao_recon_monthly[0, y] = np.nansum(sst_w[y]*gto[0])/20  # 2beta
            nao_recon_monthly[1, y] = np.nansum(sst_s[y]*gto[1])/20

        if months == 1:
            for y in range(12):
                nao_recon[0, y] = nao_recon_monthly[0, months*y]
                nao_recon[1, y] = nao_recon_monthly[1, months*y]
        else:
            for y in range(12):
                nao_recon[0, y] = np.mean(nao_recon_monthly[0, months*y:months*(y+1)])
                nao_recon[1, y] = np.mean(nao_recon_monthly[1, months*y:months*(y+1)])
        nao_recons[s] = nao_recon
    mn = np.mean(nao_recons[0], axis=1)
    sd = np.std(nao_recons[0], axis=1)
    for i in range(4):
        for j in range(2):
            nao_recons[i, j] = (nao_recons[i, j]-mn[j])/sd[j]
    return nao_recons


def sst_means():
    lat = np.linspace(89.5, -89.5, 180)
    lon = np.linspace(-179.5, 179.5, 360)
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc', 'sst')[:1776]
    sst_glob_month = np.zeros((1776))
    sst_amo_month = np.zeros((1776))
    sst_ind_month = np.zeros((1776))
    for i in range(1776):
        test = np.where(sst[i] < -999, 0, 1)
        sst1 = np.where(sst[i] < -999, 0, sst[i])
        sst_glob_month[i] = np.sum((sst1*meshlatweight)[30:150])/(np.sum((meshlatweight*test)[30:150]))
        sst_amo_month[i] = np.sum((sst1*meshlatweight)[30:90, 100:180])/(np.sum((meshlatweight*test)[30:90, 100:180]))
        sst_ind_month[i] = np.sum((sst1*meshlatweight)[70:110, 250:280])/(np.sum((meshlatweight*test)[70:110, 250:280]))

    sst_glob_ts = np.mean(np.reshape(sst_glob_month, (-1, 12)), axis=1)
    sst_amo_ts = np.mean(np.reshape(sst_amo_month, (-1, 12)), axis=1)
    sst_ind_ts = np.mean(np.reshape(sst_ind_month, (-1, 12)), axis=1)
    sst_glob = sst_glob_ts - sst_glob_ts[31:101].mean()
    sst_amo = sst_amo_ts - sst_amo_ts[31:101].mean()
    sst_ind = sst_ind_ts - sst_ind_ts[31:101].mean()

    sst_amo_wr = sst_amo - sst_glob
    sst_ind_wr = sst_ind - sst_glob

    sst_amo_wrsl = lanczos(sst_amo_wr, wi=11, co=1/16)
    sst_ind_wrsl = lanczos(sst_ind_wr, wi=11, co=1/16)
    
    ind_trend = np.zeros((108))
    for i in range(108):
        ind_trend[i] = np.polyfit(np.arange(0, 30), sst_ind_wrsl[i:i+30], 1)[0]
    print(np.corrcoef(ind_trend,sst_amo_wrsl[15:-15])[0, 1])

##############################
#### Circulation analysis ####
##############################


def regress_circu850():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat1)])
    meshlat[:, :] = lat1
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    both = compare('item15201_daily_mean')
    lst = dbase()
    anom = np.zeros((145, 192))
    u_response = np.zeros((2, len(both), 144, 192))
    sst_index = np.zeros(len(both))
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    for i, item in enumerate(both):
        anom = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        u_djfc = ucontrol[330:420]
        u_jjac = ucontrol[510:600]
        u_djf = u[330:420]
        u_jja = u[510:600]
        u_response[0, i] = np.mean(u_djf, axis=0) - np.mean(u_djfc, axis=0)
        u_response[1, i] = np.mean(u_jja, axis=0) - np.mean(u_jjac, axis=0)
        anom = anom[::-1, :]
        test = np.where(lsm > 999, 0, 1)
        sst1 = np.where(lsm > 999, 0, anom)
        sst_index[i] = np.sum((sst1*meshlatweight)[60:97, 16:57])/(np.sum((meshlatweight*test)[60:97, 16:57]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    regmap_smoothed = np.zeros((2, 144, 192))
    regmap_smoothed_sig = np.zeros((2, 144, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(144):
        for j in range(192):
            for a in range(2):
                regmap_smoothed[a, i, j] = np.cov(sst_index, u_response[a, :, i, j])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * np.cos(lat1[i]) * (np.pi/180)**2)
                r = np.cov(sst_index, u_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(sst_index, u_response[a, :, i, j])[1, 1]*np.cov(sst_index, u_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i, j] = regmap_smoothed[a, i, j]
    return regmap_smoothed, regmap_smoothed_sig


def regress_circvars():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat1)])
    meshlat[:, :] = lat1
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    both = compare('item15201_monthly_mean')
    lst = dbase()
    anom = np.zeros((145, 192))
    u_response = np.zeros((2, len(both), 144, 192))
    v_response = np.zeros((2, len(both), 144, 192))
    p_response = np.zeros((2, len(both), 145, 192))
    ssti_index = np.zeros(len(both))
    sstp_index = np.zeros(len(both))
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    for i, item in enumerate(both):
        anom = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        vcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        v = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        pcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        p = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        u_djfc = ucontrol[11:14]
        u_jjac = ucontrol[17:20]
        u_djf = u[11:14]
        u_jja = u[17:20]
        u_response[0, i] = np.mean(u_djf, axis=0) - np.mean(u_djfc, axis=0)
        u_response[1, i] = np.mean(u_jja, axis=0) - np.mean(u_jjac, axis=0)
        v_djfc = vcontrol[11:14]
        v_jjac = vcontrol[17:20]
        v_djf = v[11:14]
        v_jja = v[17:20]
        v_response[0, i] = np.mean(v_djf, axis=0) - np.mean(v_djfc, axis=0)
        v_response[1, i] = np.mean(v_jja, axis=0) - np.mean(v_jjac, axis=0)
        p_djfc = pcontrol[11:14]
        p_jjac = pcontrol[17:20]
        p_djf = p[11:14]
        p_jja = p[17:20]
        p_response[0, i] = np.mean(p_djf, axis=0) - np.mean(p_djfc, axis=0)
        p_response[1, i] = np.mean(p_jja, axis=0) - np.mean(p_jjac, axis=0)
        anom = anom[::-1, :]
        test = np.where(lsm > 999, 0, 1)
        sst1 = np.where(lsm > 999, 0, anom)
        #large ind
        #sst_index[i] = np.sum((sst1*meshlatweight)[60:97, 16:57])/(np.sum((meshlatweight*test)[60:97, 16:57]))
        # small ind
        # ssti_index[i] = np.sum((sst1*meshlatweight)[72:89, 24:49])/(np.sum((meshlatweight*test)[72:89, 24:49]))
        # ATL
        ssti_index[i] = np.sum((sst1*meshlatweight)[60:77, 160:])/(np.sum((meshlatweight*test)[60:77, 160:]))
        # pac
        # sstp_index[i] = np.sum((sst1*meshlatweight)[72:89, 80:105])/(np.sum((meshlatweight*test)[72:89, 80:105]))
        # north pac
        sstp_index[i] = np.sum((sst1*meshlatweight)[52:65, 80:105])/(np.sum((meshlatweight*test)[52:65, 80:105]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    gtoui = np.zeros((2, 144, 192))
    gtoui_sig = np.zeros((2, 144, 192))
    gtoui_sig[:] = np.NAN
    gtovi = np.zeros((2, 144, 192))
    gtovi_sig = np.zeros((2, 144, 192))
    gtovi_sig[:] = np.NAN
    gtoup = np.zeros((2, 144, 192))
    gtoup_sig = np.zeros((2, 144, 192))
    gtoup_sig[:] = np.NAN
    gtovp = np.zeros((2, 144, 192))
    gtovp_sig = np.zeros((2, 144, 192))
    gtovp_sig[:] = np.NAN
    gtopi = np.zeros((2, 145, 192))
    gtopi_sig = np.zeros((2, 145, 192))
    gtopi_sig[:] = np.NAN
    gtopp = np.zeros((2, 145, 192))
    gtopp_sig = np.zeros((2, 145, 192))
    gtopp_sig[:] = np.NAN
    for i in range(144):
        for j in range(192):
            for a in range(2):
                gtoui[a, i, j] = np.cov(ssti_index, u_response[a, :, i, j])[0, 1]/(np.cov(ssti_index, u_response[a, :, i, j])[0, 0])
                r = np.cov(ssti_index, u_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, u_response[a, :, i, j])[1, 1]*np.cov(ssti_index, u_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    gtoui_sig[a, i, j] = gtoui[a, i, j]
                gtovi[a, i, j] = np.cov(ssti_index, v_response[a, :, i, j])[0, 1]/(np.cov(ssti_index, v_response[a, :, i, j])[0, 0])
                r = np.cov(ssti_index, v_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, v_response[a, :, i, j])[1, 1]*np.cov(ssti_index, v_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    gtovi_sig[a, i, j] = gtovi[a, i, j]
                gtoup[a, i, j] = np.cov(sstp_index, u_response[a, :, i, j])[0, 1]/(np.cov(sstp_index, u_response[a, :, i, j])[0, 0])
                r = np.cov(sstp_index, u_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, u_response[a, :, i, j])[1, 1]*np.cov(sstp_index, u_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    gtoup_sig[a, i, j] = gtoup[a, i, j]
                gtovp[a, i, j] = np.cov(sstp_index, v_response[a, :, i, j])[0, 1]/(np.cov(sstp_index, v_response[a, :, i, j])[0, 0])
                r = np.cov(sstp_index, v_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, v_response[a, :, i, j])[1, 1]*np.cov(sstp_index, v_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    gtovp_sig[a, i, j] = gtovp[a, i, j]
    for i in range(145):
        for j in range(192):
            for a in range(2):
                gtopi[a, i, j] = np.cov(ssti_index, p_response[a, :, i, j])[0, 1]/(np.cov(ssti_index, p_response[a, :, i, j])[0, 0])
                r = np.cov(ssti_index, p_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, p_response[a, :, i, j])[1, 1]*np.cov(ssti_index, p_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    gtopi_sig[a, i, j] = gtopi[a, i, j]
                gtopp[a, i, j] = np.cov(sstp_index, p_response[a, :, i, j])[0, 1]/(np.cov(sstp_index, p_response[a, :, i, j])[0, 0])
                r = np.cov(sstp_index, p_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, p_response[a, :, i, j])[1, 1]*np.cov(sstp_index, p_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    gtopp_sig[a, i, j] = gtopp[a, i, j]
    return gtoui, gtoui_sig, gtovi, gtovi_sig, gtoup, gtoup_sig, gtovp, gtovp_sig, gtopi, gtopi_sig, gtopp, gtopp_sig


def regress_circvars_summer():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat1)])
    meshlat[:, :] = lat1
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    both = compare('item15201_monthly_mean')
    lst = dbase()
    anom = np.zeros((145, 192))
    u_response = np.zeros((len(both), 144, 192))
    v_response = np.zeros((len(both), 144, 192))
    p_response = np.zeros((len(both), 145, 192))
    pr_response = np.zeros((len(both), 145, 192))
    ssti_index = np.zeros(len(both))
    sstp_index = np.zeros(len(both))
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    for i, item in enumerate(both):
        anom = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        vcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        v = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        pcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        p = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        prcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item5216_monthly_mean/item5216_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item5216_monthly_mean')[:, 0, :]
        pr = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item5216_monthly_mean/item5216_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item5216_monthly_mean')[:, 0, :]
        u_jjac = ucontrol[17:20]
        u_jja = u[17:20]
        u_response[i] = np.mean(u_jja, axis=0) - np.mean(u_jjac, axis=0)
        v_jjac = vcontrol[17:20]
        v_jja = v[17:20]
        v_response[i] = np.mean(v_jja, axis=0) - np.mean(v_jjac, axis=0)
        p_jjac = pcontrol[17:20]
        p_jja = p[17:20]
        p_response[i] = np.mean(p_jja, axis=0) - np.mean(p_jjac, axis=0)
        pr_jjac = prcontrol[17:20]
        pr_jja = pr[17:20]
        pr_response[i] = np.mean(pr_jja, axis=0) - np.mean(pr_jjac, axis=0)
        anom = anom[::-1, :]
        test = np.where(lsm > 999, 0, 1)
        sst1 = np.where(lsm > 999, 0, anom)
        # pac
        #ssti_index[i] = np.sum((sst1*meshlatweight)[48:65, 54:81])/(np.sum((meshlatweight*test)[48:65, 54:81]))
        # pac dipole
        #sstp_index[i] = np.sum((sst1*meshlatweight)[48:65, 64:110])/(np.sum((meshlatweight*test)[48:65, 64:110])) - np.sum((sst1*meshlatweight)[65:85, 64:110])/(np.sum((meshlatweight*test)[65:85, 64:110]))
        # ind
        ssti_index[i] = np.sum((sst1*meshlatweight)[48:69, 24:54])/(np.sum((meshlatweight*test)[48:69, 24:54]))
        # pac neg
        sstp_index[i] = np.sum((sst1*meshlatweight)[65:85, 75:118])/(np.sum((meshlatweight*test)[65:85, 75:118]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    gtou = np.zeros((2, 144, 192))
    gtou_sig = np.zeros((2, 144, 192))
    gtou_sig[:] = np.NAN
    gtov = np.zeros((2, 144, 192))
    gtov_sig = np.zeros((2, 144, 192))
    gtov_sig[:] = np.NAN
    gtop = np.zeros((2, 145, 192))
    gtop_sig = np.zeros((2, 145, 192))
    gtop_sig[:] = np.NAN
    gtopr = np.zeros((2, 145, 192))
    gtopr_sig = np.zeros((2, 145, 192))
    gtopr_sig[:] = np.NAN
    for i in range(144):
        for j in range(192):
            gtou[0, i, j] = np.cov(ssti_index, u_response[:, i, j])[0, 1]/(np.cov(ssti_index, u_response[:, i, j])[0, 0])
            r = np.cov(ssti_index, u_response[:, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, u_response[:, i, j])[1, 1]*np.cov(ssti_index, u_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtou_sig[0, i, j] = gtou[0, i, j]
            gtov[0, i, j] = np.cov(ssti_index, v_response[:, i, j])[0, 1]/(np.cov(ssti_index, v_response[:, i, j])[0, 0])
            r = np.cov(ssti_index, v_response[:, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, v_response[:, i, j])[1, 1]*np.cov(ssti_index, v_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtov_sig[0, i, j] = gtov[0, i, j]
            gtou[1, i, j] = np.cov(sstp_index, u_response[:, i, j])[0, 1]/(np.cov(sstp_index, u_response[:, i, j])[0, 0])
            r = np.cov(sstp_index, u_response[:, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, u_response[:, i, j])[1, 1]*np.cov(sstp_index, u_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtou_sig[1, i, j] = gtou[1, i, j]
            gtov[1, i, j] = np.cov(sstp_index, v_response[:, i, j])[0, 1]/(np.cov(sstp_index, v_response[:, i, j])[0, 0])
            r = np.cov(sstp_index, v_response[:, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, v_response[:, i, j])[1, 1]*np.cov(sstp_index, v_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtov_sig[1, i, j] = gtov[1, i, j]
    for i in range(145):
        for j in range(192):
            gtop[0, i, j] = np.cov(ssti_index, p_response[:, i, j])[0, 1]/(np.cov(ssti_index, p_response[:, i, j])[0, 0])
            r = np.cov(ssti_index, p_response[:, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, p_response[:, i, j])[1, 1]*np.cov(ssti_index, p_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtop_sig[0, i, j] = gtop[0, i, j]
            gtop[1, i, j] = np.cov(sstp_index, p_response[:, i, j])[0, 1]/(np.cov(sstp_index, p_response[:, i, j])[0, 0])
            r = np.cov(sstp_index, p_response[:, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, p_response[:, i, j])[1, 1]*np.cov(sstp_index, p_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtop_sig[0, i, j] = gtop[0, i, j]
            gtopr[0, i, j] = np.cov(ssti_index, pr_response[:, i, j])[0, 1]/(np.cov(ssti_index, pr_response[:, i, j])[0, 0])
            r = np.cov(ssti_index, pr_response[:, i, j])[0, 1]/np.sqrt(np.cov(ssti_index, pr_response[:, i, j])[1, 1]*np.cov(ssti_index, pr_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtopr_sig[0, i, j] = gtopr[0, i, j]
            gtopr[1, i, j] = np.cov(sstp_index, pr_response[:, i, j])[0, 1]/(np.cov(sstp_index, pr_response[:, i, j])[0, 0])
            r = np.cov(sstp_index, pr_response[:, i, j])[0, 1]/np.sqrt(np.cov(sstp_index, pr_response[:, i, j])[1, 1]*np.cov(sstp_index, pr_response[:, i, j])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                gtopr_sig[0, i, j] = gtopr[0, i, j]
    return gtou, gtou_sig, gtov, gtov_sig, gtop, gtop_sig, gtopr, gtopr_sig


'''
def regress_circv200():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat1)])
    meshlat[:, :] = lat1
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    both = compare('item15202_monthly_mean')
    lst = dbase()
    anom = np.zeros((145, 192))
    u_response = np.zeros((2, len(both), 144, 192))
    sst_index = np.zeros(len(both))
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    for i, item in enumerate(both):
        anom = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        u_djfc = ucontrol[11:14]
        u_jjac = ucontrol[17:20]
        u_djf = u[11:14]
        u_jja = u[17:20]
        u_response[0, i] = np.mean(u_djf, axis=0) - np.mean(u_djfc, axis=0)
        u_response[1, i] = np.mean(u_jja, axis=0) - np.mean(u_jjac, axis=0)
        anom = anom[::-1, :]
        test = np.where(lsm > 999, 0, 1)
        sst1 = np.where(lsm > 999, 0, anom)
        #large ind
        #sst_index[i] = np.sum((sst1*meshlatweight)[60:97, 16:57])/(np.sum((meshlatweight*test)[60:97, 16:57]))
        # small ind
        # sst_index[i] = np.sum((sst1*meshlatweight)[72:89, 24:44])/(np.sum((meshlatweight*test)[72:89, 24:44]))
        # pac
        sst_index[i] = np.sum((sst1*meshlatweight)[72:81, 80:97])/(np.sum((meshlatweight*test)[72:81, 80:97]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    regmap_smoothed = np.zeros((2, 144, 192))
    regmap_smoothed_sig = np.zeros((2, 144, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(144):
        for j in range(192):
            for a in range(2):
                regmap_smoothed[a, i, j] = np.cov(sst_index, u_response[a, :, i, j])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * np.cos(lat1[i]) * (np.pi/180)**2)
                r = np.cov(sst_index, u_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(sst_index, u_response[a, :, i, j])[1, 1]*np.cov(sst_index, u_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i, j] = regmap_smoothed[a, i, j]
    return regmap_smoothed, regmap_smoothed_sig


def regress_circmslp():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_b01f_1999-01_2000-12.nc', 'latitude1')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_b01f_1999-01_2000-12.nc', 'longitude1')
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat1)])
    meshlat[:, :] = lat1
    meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((145, 192))
    u_response = np.zeros((2, len(both), 145, 192))
    sst_index = np.zeros(len(both))
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    for i, item in enumerate(both):
        anom = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        u_djfc = ucontrol[11:14]
        u_jjac = ucontrol[17:20]
        u_djf = u[11:14]
        u_jja = u[17:20]
        u_response[0, i] = np.mean(u_djf, axis=0) - np.mean(u_djfc, axis=0)
        u_response[1, i] = np.mean(u_jja, axis=0) - np.mean(u_jjac, axis=0)
        anom = anom[::-1, :]
        test = np.where(lsm > 999, 0, 1)
        sst1 = np.where(lsm > 999, 0, anom)
        #large ind
        #sst_index[i] = np.sum((sst1*meshlatweight)[60:97, 16:57])/(np.sum((meshlatweight*test)[60:97, 16:57]))
        # small ind
        # sst_index[i] = np.sum((sst1*meshlatweight)[72:89, 24:44])/(np.sum((meshlatweight*test)[72:89, 24:44]))
        # pac
        sst_index[i] = np.sum((sst1*meshlatweight)[72:81, 80:97])/(np.sum((meshlatweight*test)[72:81, 80:97]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    regmap_smoothed = np.zeros((2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(145):
        for j in range(192):
            for a in range(2):
                regmap_smoothed[a, i, j] = np.cov(sst_index, u_response[a, :, i, j])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * np.cos(lat1[i]) * (np.pi/180)**2)
                r = np.cov(sst_index, u_response[a, :, i, j])[0, 1]/np.sqrt(np.cov(sst_index, u_response[a, :, i, j])[1, 1]*np.cov(sst_index, u_response[a, :, i, j])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i, j] = regmap_smoothed[a, i, j]
    return regmap_smoothed, regmap_smoothed_sig
'''


def clim_rpm():
    from scipy import stats
    both = compare('item15202_monthly_mean')
    lst = dbase()
    u_c = np.zeros((2, 144, 192))
    v_c = np.zeros((2, 144, 192))
    for i, item in enumerate(both):
        vcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_monthly_mean/item15201_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_monthly_mean')[:, 1, :]
        u_c[0] += np.mean(ucontrol[11:14], axis=0)
        u_c[1] += np.mean(ucontrol[17:20], axis=0)
        v_c[0] += np.mean(vcontrol[11:14], axis=0)
        v_c[1] += np.mean(vcontrol[17:20], axis=0)
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    u_c /= 5544
    v_c /= 5544
    return u_c, v_c


def clim_rpm_u850():
    import glob
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/*')
    u_c = np.zeros((2, 144, 192))
    j = 0
    for i in a:
        ucontrol = ncread(i, 'item15201_daily_mean')[:, 0, :]
        u_c[0] += np.mean(ucontrol[330:420], axis=0)
        u_c[1] += np.mean(ucontrol[510:630], axis=0)
        j += 1
        print('Done: ' + str(j) + ' out of ' + str(len(a)))
    u_c /= len(a)
    return u_c

'''
Think this is incorrect as don't want to use non control data!
def deser_circv200():
    import scipy.stats as stats
    def nao_indivyear(eof, indiv_years, lat, lon, ms='no'):
        # project responses onto eof
        proj = eof_response(indiv_years, eof, lat, lon)
        proj_m = np.mean(proj)
        proj_s = np.std(proj)
        proj = (proj-proj_m)/proj_s
        if ms == 'yes':
            return proj, proj_m, proj_s
        return proj
    n_eof = 1
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'longitude0')
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    both = compare('item15202_monthly_mean')
    mslp_w = np.zeros((len(both), 145, 192))
    mslp_s = np.zeros((len(both), 145, 192))
    mslp_wf = np.zeros((len(both), 145, 192))
    mslp_sf = np.zeros((len(both), 145, 192))
    u_response = np.zeros((2, len(both), 144, 192))
    for i, item in enumerate(both):
        mslp_96_c = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :] / 100
        mslp_96_f = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :] / 100
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        u_djfc = ucontrol[11:14]
        u_jjac = ucontrol[17:20]
        u_djf = u[11:14]
        u_jja = u[17:20]
        u_response[0, i] = np.mean(u_djf, axis=0) - np.mean(u_djfc, axis=0)
        u_response[1, i] = np.mean(u_jja, axis=0) - np.mean(u_jjac, axis=0)
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        field1 = mslp_96_c[:, :, 96:]
        field2 = mslp_96_c[:, :, :96]
        field3 = mslp_96_f[:, :, 96:]
        field4 = mslp_96_f[:, :, :96]
        mslp_c = np.concatenate((field1, field2), axis=2)
        mslp_f = np.concatenate((field3, field4), axis=2)
        mslp_w[i] = np.mean(mslp_c[windices], axis=0)
        mslp_s[i] = np.mean(mslp_c[sindices], axis=0)
        mslp_wf[i] = np.mean(mslp_f[windices], axis=0)
        mslp_sf[i] = np.mean(mslp_f[sindices], axis=0)
        print(str(i))
    eofs = np.zeros((2, 145, 192))
    var = np.zeros((2, 3))
    nao = np.zeros((2, len(mslp_w)))
    nao_f = np.zeros((2, len(mslp_w)))
    eofs[0], var[0] = eof_clim(mslp_w, lat_n96, lon_n96, n_eof)
    nao[0], ms, sd = nao_indivyear(eofs[0, 8:57, 48:118], mslp_w[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118], ms='yes')
    proj = eof_response(mslp_wf[:, 8:57, 48:118], eofs[0, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    nao_f[0] = (proj-ms)/sd
    eofs[1], var[1] = eof_clim(mslp_s, lat_n96, lon_n96, n_eof)
    nao[1], ms, sd = nao_indivyear(eofs[1, 8:57, 48:118], mslp_s[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118], ms='yes')
    proj = eof_response(mslp_sf[:, 8:57, 48:118], eofs[1, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    nao_f[1] = (proj-ms)/sd

    nao_r = nao_f - nao

    v200_nao = np.zeros((2, 144, 192))
    for i in range(144):
        for j in range(192):
            for a in range(2):
                v200_nao[a, i, j] = stats.linregress(nao_r[a], u_response[a, :, i, j])[0]
    return v200_nao
'''


def deser_circv200():
    import scipy.stats as stats

    def nao_indivyear(eof, indiv_years, lat, lon, ms='no'):
        # project responses onto eof
        proj = eof_response(indiv_years, eof, lat, lon)
        proj_m = np.mean(proj)
        proj_s = np.std(proj)
        proj = (proj-proj_m)/proj_s
        if ms == 'yes':
            return proj, proj_m, proj_s
        return proj

    n_eof = 1
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    both = compare('item15202_monthly_mean')
    mslp_w = np.zeros((len(both), 145, 192))
    mslp_s = np.zeros((len(both), 145, 192))
    u_response = np.zeros((2, len(both), 144, 192))
    for i, item in enumerate(both):
        mslp_96_c = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :] / 100
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        u_response[0, i] = ucontrol[11:14].mean(axis=0)
        u_response[1, i] = ucontrol[17:20].mean(axis=0)
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        field1 = mslp_96_c[:, :, 96:]
        field2 = mslp_96_c[:, :, :96]
        mslp_c = np.concatenate((field1, field2), axis=2)
        mslp_w[i] = np.mean(mslp_c[windices], axis=0)
        mslp_s[i] = np.mean(mslp_c[sindices], axis=0)
        print(str(i))
    eofs = np.zeros((2, 145, 192))
    var = np.zeros((2, 3))
    nao = np.zeros((2, len(mslp_w)))
    eofs[0], var[0] = eof_clim(mslp_w, lat_n96, lon_n96, n_eof)
    nao[0] = nao_indivyear(eofs[0, 8:57, 48:118], mslp_w[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    eofs[1], var[1] = eof_clim(mslp_s, lat_n96, lon_n96, n_eof)
    nao[1] = nao_indivyear(eofs[1, 8:57, 48:118], mslp_s[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])

    v200_nao = np.zeros((2, 144, 192))
    for i in range(144):
        for j in range(192):
            for a in range(2):
                v200_nao[a, i, j] = stats.linregress(nao[a], u_response[a, :, i, j])[0]
    return v200_nao


#####################
####### CESM ########
#####################


def cesm_eof(n_eof=1):
    from scipy import interpolate
    import glob
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.arange(0, 360, 360/192)
    a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/RPM/CESM-LE/coupled_control/psl/*'))
    mslp = np.zeros((1200*18+12, 192, 288))
    for i in range(len(a)):
        mslp[i*1200:(i+1)*1200] = ncread(a[i], 'PSL')[:1200] / 100
    mslp[-12:] = ncread(a[-1], 'PSL')[1200:] / 100
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/CESM-LE/coupled_control/psl/b.e11.B1850C5CN.f09_g16.005.cam.h0.PSL.040001-049912.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/RPM/CESM-LE/coupled_control/psl/b.e11.B1850C5CN.f09_g16.005.cam.h0.PSL.040001-049912.nc', 'lon')

    mslp_96 = np.zeros((np.shape(mslp)[0], 145, 192))
    for i in range(np.shape(mslp)[0]):
        f = interpolate.interp2d(lon, lat[::-1], mslp[i])
        mslp_96[i] = f(lon_n96, lat_n96)

    sindices = np.zeros((3*(np.shape(mslp_96)[0]//12-1)))
    windices = np.zeros((3*(np.shape(mslp_96)[0]//12-1)))
    for i in range((np.shape(mslp_96)[0]//12-1)):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    field1 = mslp_96[:, :, 96:]
    field2 = mslp_96[:, :, :96]
    mslp_c = np.concatenate((field1, field2), axis=2)
    mslp_w = mslp_c[windices]
    mslp_s = mslp_c[sindices]
    mslp_w = np.add.reduceat(mslp_w, range(0, len(mslp_w), 3))/3
    mslp_s = np.add.reduceat(mslp_s, range(0, len(mslp_s), 3))/3

    eofs = np.zeros((2, 145, 192))
    var = np.zeros((2, 3))
    nao = np.zeros((2, len(mslp_w)))
    eofs[0], var[0] = eof_clim(mslp_w, lat_n96, lon_n96, n_eof)
    nao[0] = nao_indivyear(eofs[0, 8:57, 48:118]*-1, mslp_w[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    eofs[1], var[1] = eof_clim(mslp_s, lat_n96, lon_n96, n_eof)
    nao[1] = nao_indivyear(eofs[1, 8:57, 48:118], mslp_s[:, 8:57, 48:118], lat_n96[8:57], lon_n96[48:118])
    return eofs, var, nao


def reconstruct_nao_cesm(gto, sst, region='glob'):
    from scipy import interpolate
    from mpl_toolkits.basemap import shiftgrid
    lat_n96 = np.linspace(90, -90, 145)
    lon_n96 = np.linspace(0, 358.125, 192)
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    lon_sst_shift = np.linspace(.5, 359.5, 360)

    sindices = np.zeros((1800*3))
    windices = np.zeros((1800*3))
    for i in range(1800):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    gto_interp = np.zeros((2, 180, 360))

    sst_0 = np.where(sst > 999, 0, sst)
    gto_0 = np.where(np.isnan(gto), 0, gto)
    if region == 'atl':
        gto_0[:, :, 16:145] = 0 #90-0 for basin
    elif region == 'ind':
        gto_0[:, :, :16] = 0
        gto_0[:, :, 59:] = 0
    elif region == 'pac':
        gto_0[:, :, :59] = 0 #no equals fo basin
        gto_0[:, :, 145:] = 0 # 240E for basin
    for i in range(2):
        f = interpolate.interp2d(lon_n96, lat_n96[::-1], gto_0[i, :])
        gto_interp[i, :] = f(lon_sst_shift, lat_sst)
        gto_interp[i], lon1 = shiftgrid(181., gto_interp[i], lon_sst_shift,
                                        start=False)

    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2
    nao_recon_monthly = np.zeros((2, 3*1800))
    nao_recon = np.zeros((2, 1800))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight

    for y in range(1800*3):
        nao_recon_monthly[0, y] = np.nansum(sst_w[y]*gto_interp[0])/20  # 2beta
        nao_recon_monthly[1, y] = np.nansum((sst_s[y])*gto_interp[1])/20

    for y in range(1800):
        nao_recon[0, y] = np.mean(nao_recon_monthly[0, 3*y:3*(y+1)])
        nao_recon[1, y] = np.mean(nao_recon_monthly[1, 3*y:3*(y+1)])
    return nao_recon


def n34ind(sst):
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst = np.linspace(-179.5, 179.5, 360)
    sindices = np.zeros((1800*3))
    for i in range(1800):
        sindices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)

    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180)

    sst_0 = np.where(sst > 999, 0, sst)
    mask = np.copy(sst)[0]
    mask = np.where(mask > 999, 0, 1)
    sst_1 = sst_0[sindices]
    sst_1 = np.mean(np.reshape(sst_1, (-1, 3, 180, 360)), axis=1)
    n34 = np.mean((sst_1*meshlatweight)[:, 85:95, 10:60], axis=(1, 2)) / np.mean((meshlatweight*mask)[85:95, 10:60])
    n34 = (n34 - n34.mean())
    return n34


def load_ssts():
    import glob
    sst = np.zeros((1200*18+12, 180, 360))
    a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/RPM/CESM-LE/coupled_control/tos/rotated/*'))
    for i in range(len(a)):
        sst[i*1200:(i+1)*1200] = ncread(a[i], 'SST')[:1200, 0]
        print(str(i))
    sst[-12:] = ncread(a[-1], 'SST')[1200:, 0]
    return sst


################################
### Pacific stuff for Chris ####
################################


def regress_pac():
    from scipy import stats
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(-180, 180, 1.875)
    lat1 = np.linspace(90, -90, 145) * np.pi/180
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    jec = np.zeros((len(both)))
    je = np.zeros((len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]/100
        u_djfc = ucontrol[11:14].mean(axis=0)
        u_djf = u[11:14].mean(axis=0)
        jec[i] = (np.mean((u_djfc*meshlatweight)[20:49, 96:129]) / np.mean(meshlatweight[20:49, 96:129]))
        je[i] = (np.mean((u_djf*meshlatweight)[20:49, 96:129]) / np.mean(meshlatweight[20:49, 96:129]))
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))

    anom = anom[:, ::-1, :]
    je_d = je - jec
    regmap_smoothed = np.zeros((145, 192))
    regmap_smoothed_sig = np.zeros((145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            regmap_smoothed[i+23, j] = np.cov(anom[:, i+23, j], je_d[:])[0, 1]*3/(4 * 6371**2 * 1.25 * 1.875 * np.cos(lat1[i+23]) * (np.pi/180)**2)
            r = np.cov(anom[:, i+23, j], je_d[:])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], je_d[:])[1, 1]*np.cov(anom[:, i+23, j], je_d[:])[0, 0])
            t = r*np.sqrt((len(both)-2)/(1-r**2))
            p = 1 - stats.norm.cdf(np.abs(t))
            sig = np.greater_equal(5, p*100*2).astype(int)
            if sig == 1:
                regmap_smoothed_sig[i+23, j] = regmap_smoothed[i+23, j]
    return regmap_smoothed, regmap_smoothed_sig
