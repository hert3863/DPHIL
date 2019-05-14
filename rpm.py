# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:27:20 2017

@author: bakerh

USE LAT 90 to MINUS 90 WHEN PLOTTING
"""

import numpy as np


def ncread(filelocation, invariable):
    '''
    ncread outputs numpy arrays of called variables in netCDF4 file

    Parameters
    ----------
    fileloaction : str
        NetCDFfile
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variable : array
        Array of variable specified
    '''
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    nc_f = filelocation  # filename
    nc_fid = Dataset(nc_f, 'r')
    variable = nc_fid.variables[invariable][:]
    return variable


def shiftlat(grid):
    latold = np.arange(89.375, -90.625, -1.25)
    latnew = np.arange(90, -91.25, -1.25)
    regrid = np.zeros((np.ma.size(grid, axis=0), 145, 192))
    for i in range(143):
        regrid[:, i+1, :] = ((grid[:, i, :]*np.cos(latold[i]*np.pi/180) +
                              grid[:, i+1, :]*np.cos(latold[i+1]*np.pi/180)) /
                             (2*np.cos(latnew[i+1]*np.pi/180)))
    return regrid


def jetdist(u850all, latall):
    """
    Determines the jet distributions

    Parameters
    ----------
    u850all: array
        wind data at u850
    lat: array
        latitudes of data

    Returns
    -------
    wdist: array
        lat distribution
    wweights: array
        speeds at each lat
    """
    u850 = np.mean(u850all[15:61, 160:], axis=1)
    lat = latall[15:61]
    wwmx = np.max(u850)
    windex = np.argmax(u850)
    wdmax = lat[windex]
    wdist = [wdmax]
    wweights = [wwmx]

    windex -= 1
    while u850[windex+1] >= u850[windex] and u850[windex+1] >= 0 and windex >= 0:
        if u850[windex] < 0:
            lat_0 = lat[windex] - u850[windex] * ((lat[windex+1]-lat[windex]) /
                                                  (u850[windex+1] -
                                                   u850[windex]))
            wdist.insert(0, lat_0)
            wweights.insert(0, 0)
            windex -= 1
        else:
            wdist.insert(0, lat[windex])
            wweights.insert(0, u850[windex])
            windex -= 1

    windex = np.argmax(u850)
    windex += 1
    while u850[windex-1] >= u850[windex] and u850[windex-1] >= 0:
        if u850[windex] < 0:
            lat_0 = lat[windex-1] - u850[windex-1] * ((lat[windex] -
                                                       lat[windex-1]) /
                                                      (u850[windex] -
                                                       u850[windex-1]))
            wdist.append(lat_0)
            wweights.append(0)
            windex += 1
        else:
            wdist.append(lat[windex])
            wweights.append(u850[windex])
            windex += 1

    return wdist, wweights


def jetind_poly(u850, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    HadAM3P output

    Parameters
    ----------
    u850: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    jetspeed: float
        max jet speed
    jetlatitude: float
        position of maximum
    """
    # from sympy import Symbol
    # from sympy import solve
    from scipy import optimize
    speeds = np.zeros((np.ma.size(u850, axis=0)))
    lats = np.zeros((np.ma.size(u850, axis=0)))
    for t in range(np.ma.size(u850, axis=0)):
        wd, ww = jetdist(u850[t, :], lat)
        wm = np.argmax(ww)
        wp = np.polyfit(wd[wm-1:wm+2], ww[wm-1:wm+2], 2)

        def f(x): return wp[0] * x**2 + wp[1] * x + wp[2]
        xmax = optimize.fmin(lambda x: -f(x), wm, disp=False)
        ymax = f(xmax)

        speeds[t] = ymax
        lats[t] = xmax[0]

    return speeds, lats


def jetind(u850, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    HadAM3P output

    Parameters
    ----------
    u850: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    jetspeed: float
        max jet speed
    jetlatitude: float
        position of maximum
    """
    # from sympy import Symbol
    # from sympy import solve
    speeds = np.zeros((np.ma.size(u850, axis=0)))
    lats = np.zeros((np.ma.size(u850, axis=0)))
    u850_mean = np.mean(u850[:, 15:60, 160:], axis=2)
    lat_dist = lat[15:60]
    for t in range(np.ma.size(u850, axis=0)):
        speeds[t] = np.max(u850_mean[t, :])
        lats[t] = lat_dist[np.argmax(u850_mean[t, :])]

    return speeds, lats


def jetindR(u850, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    HadAM3P output

    Parameters
    ----------
    u850: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    jetspeed: float
        max jet speed
    jetlatitude: float
        position of maximum
    """
    # from sympy import Symbol
    # from sympy import solve
    speeds = np.zeros((np.ma.size(u850, axis=0)))
    lats = np.zeros((np.ma.size(u850, axis=0)))
    u850_mean = np.mean(u850[:, 7:31, 120:], axis=2)
    lat_dist = lat[7:31]
    for t in range(np.ma.size(u850, axis=0)):
        speeds[t] = np.max(u850_mean[t, :])
        lats[t] = lat_dist[np.argmax(u850_mean[t, :])]

    return speeds, lats


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


def ttest_n(series1, series2, siglevel=5, testtype='two'):
    """
    Determines the significance of the applied forcing

    Parameters
    ----------
    series1: array
        control run
    series2: array
        forced run

    Returns
    -------
    sig: float
        significance of difference
    """
    import scipy.stats as st
    if testtype == 'two':
        a = 2
    elif testtype == 'one':
        a = 1
    var1 = np.var(series1, ddof=1)
    var2 = np.var(series2, ddof=1)
    rho1 = np.corrcoef(series1[:-1], series1[1:])
    rho1 = abs(rho1[1, 0])
    rho2 = np.corrcoef(series2[:-1], series2[1:])
    rho2 = abs(rho2[1, 0])
    neff1 = np.ma.size(series1) * (1-rho1) / (1+rho1)
    neff2 = np.ma.size(series2) * (1-rho2) / (1+rho2)
    z = ((np.mean(series1) - np.mean(series2)) /
         (np.sqrt((var1 / neff1) + (var2 / neff2))))
    p = 1 - st.norm.cdf(np.abs(z))
    sig = np.greater_equal(siglevel, p*100*a).astype(int)

    return sig


def regress():
    from scipy import stats
    lat = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_b01f_1999-01_2000-12.nc', 'latitude0')
    lat1 = np.arange(90, -91.25, -1.25) * np.pi/180
    both = compare('item15201_daily_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    jetinds = np.zeros((2, 2, len(both)))
    jetindscontrol = np.zeros((2, 2, len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        ucontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        u = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15201_daily_mean/item15201_daily_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15201_daily_mean')[:, 0, :]
        jetinds_daily = np.zeros((2, 720))
        jetinds_dailycontrol = np.zeros((2, 720))
        jetinds_daily[0, :], jetinds_daily[1, :] = jetind(u, lat)
        jetinds_dailycontrol[0, :], jetinds_dailycontrol[1, :] = jetind(ucontrol, lat)
        jetinds[0, :, i] = np.mean(np.concatenate((jetinds_daily[:, :60],
                                   jetinds_daily[:, 330:420],
                                   jetinds_daily[:, 690:]), axis=1), axis=1)
        jetinds[1, :, i] = np.mean(np.concatenate((jetinds_daily[:, 150:240],
                                   jetinds_daily[:, 510:600]), axis=1), axis=1)
        jetindscontrol[0, :, i] = np.mean(np.concatenate((jetinds_dailycontrol[:, :60],
                                          jetinds_dailycontrol[:, 330:420],
                                          jetinds_dailycontrol[:, 690:]),
                                          axis=1), axis=1)
        jetindscontrol[1, :, i] = np.mean(np.concatenate((jetinds_dailycontrol[:, 150:240],
                                          jetinds_dailycontrol[:, 510:600]),
                                          axis=1), axis=1)
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    
    anom = anom[:, ::-1 ,:]
    
    # meshlat = np.zeros([np.ma.size(lat1), 192])
    # meshlat[:, :] = np.expand_dims(lat1, axis=1)
    # meshlatweight = np.cos(meshlat * np.pi/180) * 6371000**2 * 1.25 * 1.875
    # anom = anom * meshlatweight
    jetind_d = jetinds - jetindscontrol
    regmap = np.zeros((2, 2, 145, 192))
    regmap_p = np.ones((2, 2, 145, 192))
    regmap_sig = np.zeros((2, 2, 145, 192))
    regmap_sig[:] = np.NAN

    regmap_smoothed = np.zeros((2, 2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for a in range(2):
                for b in range(2):
                    # regmap[a, b, i+23, j], regmap_p[a, b, i+23, j] = stats.linregress(anom[:, i+23, j], jetind_d[a, b, :])[::3]
                    # if regmap_p[a, b, i+23, j] <= 0.025:  # change this if go back to lin reg sig
                    #     regmap_sig[a, b, i+23, j] = regmap[a, b, i+23, j]
                    regmap_smoothed[a, b, i+23, j] = np.cov(anom[:, i+23, j], jetind_d[a, b, :])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * np.cos(lat1[i+23]) * (np.pi/180)**2)
                    r = np.cov(anom[:, i+23, j], jetind_d[a, b, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], jetind_d[a, b, :])[1, 1]*np.cov(anom[:, i+23, j], jetind_d[a, b, :])[0, 0])
                    t = r*np.sqrt((len(both)-2)/(1-r**2))
                    p = 1 - stats.norm.cdf(np.abs(t))
                    sig = np.greater_equal(5, p*100*2).astype(int)
                    if sig == 1:
                        regmap_smoothed_sig[a, b, i+23, j] = regmap_smoothed[a, b, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def regressmode(mode='NAO'):
    # p = 0.25*(z500[:, 56, 107]-z500[:, 36, 104]+z500[:, 25, 131]-z500[:, 48, 147])
    from scipy import stats
    lat1 = np.arange(90, -91.25, -1.25)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    NAOinds = np.zeros((2, len(both)))
    slp_all = np.zeros((2, 24*len(both)))
    slp_control_all = np.zeros((2, 24*len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        slpcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        slp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        slp_control_all[:, i*24:(i+1)*24] = slpcontrol[:, 41, 187], slpcontrol[:, 20, 180]
        slp_all[:, i*24:(i+1)*24] = slp[:, 41, 187], slp[:, 20, 180]
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    slp_all = np.transpose(slp_all)
    slp_control_all = np.transpose(slp_control_all)
    for i in range(12):
        slp_all[i::12, :] = ((slp_all[i::12, :] -
                              np.mean(slp_control_all[i::12, :], axis=0)) /
                             np.std(slp_control_all[i::12, :], axis=0))
        slp_control_all[i::12, :] = ((slp_control_all[i::12, :] -
                                      np.mean(slp_control_all[i::12, :],
                                              axis=0)) /
                                     np.std(slp_control_all[i::12, :],
                                            axis=0))
    NAOind = slp_all[:, 0] - slp_all[:, 1] - (slp_control_all[:, 0] -
                                              slp_control_all[:, 1])
    for i in range(len(both)):
        NAOinds[0, i] = np.mean(np.concatenate((NAOind[24*i:24*i+2],
                                NAOind[24*i+11:24*i+14],
                                np.expand_dims(NAOind[24*i+23], axis=0)),
                                               axis=0), axis=0)
        NAOinds[1, i] = np.mean(np.concatenate((NAOind[24*i+5:24*i+8],
                                NAOind[24*i+17:24*i+20]), axis=0), axis=0)
    anom = anom[:, ::-1 ,:]
    # meshlat = np.zeros([np.ma.size(lat1), 192])
    # meshlat[:, :] = np.expand_dims(lat1, axis=1)
    # meshlatweight = np.cos(meshlat * np.pi/180) * 6371000**2 * 1.25 * 1.875
    # r_naanom = anom * meshlatweight
    regmap = np.zeros((2, 145, 192))
    regmap_p = np.ones((2, 145, 192))
    regmap_sig = np.zeros((2, 145, 192))
    regmap_sig[:] = np.NAN
    regmap_smoothed = np.zeros((2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for a in range(2):
                #regmap[a, i+23, j], regmap_p[a, i+23, j] = stats.linregress(anom[:, i+23, j], NAOinds[a, :])[::3]
                #if regmap_p[a, i+23, j] <= 0.025:
                #    regmap_sig[a, i+23, j] = regmap[a, i+23, j]
                regmap_smoothed[a, i+23, j] = np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * (np.pi/180)**2)
                r = np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], NAOinds[a, :])[1, 1]*np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i+23, j] = regmap_smoothed[a, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def regressmodeequi(mode='NAO'):
    # p = 0.25*(z500[:, 56, 107]-z500[:, 36, 104]+z500[:, 25, 131]-z500[:, 48, 147])
    from scipy import stats
    lat1 = np.arange(90, -91.25, -1.25)
    both = compare('item16222_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    NAOinds = np.zeros((2, len(both)))
    slp_all = np.zeros((2, 24*len(both)))
    slp_control_all = np.zeros((2, 24*len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        slpcontrol = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        slp = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/item16222_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item16222_monthly_mean')[:, 0, :]
        slp_control_all[:, i*24:(i+1)*24] = slpcontrol[:, 41, 187], slpcontrol[:, 20, 180]
        slp_all[:, i*24:(i+1)*24] = slp[:, 41, 187], slp[:, 20, 180]
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    slp_all = np.transpose(slp_all)
    slp_control_all = np.transpose(slp_control_all)
    for i in range(12):
        slp_all[i::12, :] = ((slp_all[i::12, :] -
                              np.mean(slp_control_all[i::12, :], axis=0)) /
                             np.std(slp_control_all[i::12, :], axis=0))
        slp_control_all[i::12, :] = ((slp_control_all[i::12, :] -
                                      np.mean(slp_control_all[i::12, :],
                                              axis=0)) /
                                     np.std(slp_control_all[i::12, :],
                                            axis=0))
    NAOind = slp_all[:, 0] - slp_all[:, 1] - (slp_control_all[:, 0] -
                                              slp_control_all[:, 1])
    for i in range(len(both)):
        NAOinds[0, i] = np.mean(np.concatenate((NAOind[24*i+2:24*i+5],
                                NAOind[24*i+14:24*i+17]), axis=0), axis=0)
        NAOinds[1, i] = np.mean(np.concatenate((NAOind[24*i+8:24*i+11],
                                NAOind[24*i+20:24*i+23]), axis=0), axis=0)
    anom = anom[:, ::-1 ,:]
    # meshlat = np.zeros([np.ma.size(lat1), 192])
    # meshlat[:, :] = np.expand_dims(lat1, axis=1)
    # meshlatweight = np.cos(meshlat * np.pi/180) * 6371000**2 * 1.25 * 1.875
    # r_naanom = anom * meshlatweight
    regmap = np.zeros((2, 145, 192))
    regmap_p = np.ones((2, 145, 192))
    regmap_sig = np.zeros((2, 145, 192))
    regmap_sig[:] = np.NAN
    regmap_smoothed = np.zeros((2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for a in range(2):
                #regmap[a, i+23, j], regmap_p[a, i+23, j] = stats.linregress(anom[:, i+23, j], NAOinds[a, :])[::3]
                #if regmap_p[a, i+23, j] <= 0.025:
                #    regmap_sig[a, i+23, j] = regmap[a, i+23, j]
                regmap_smoothed[a, i+23, j] = np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * (np.pi/180)**2)
                r = np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], NAOinds[a, :])[1, 1]*np.cov(anom[:, i+23, j], NAOinds[a, :])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i+23, j] = regmap_smoothed[a, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def regress_circ():
    from scipy import stats
    both = compare('item15202_monthly_mean')
    lst = dbase()
    anom = np.zeros((len(both), 145, 192))
    circinds = np.zeros((2, len(both)))
    wav_all = np.zeros((24*len(both)))
    wav_control_all = np.zeros((24*len(both)))
    for i, item in enumerate(both):
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[item] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        vc = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        vr = ncread('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item15202_monthly_mean/item15202_monthly_mean_' + both[i] + '_1999-01_2000-12.nc', 'item15202_monthly_mean')[:, 1, :]
        wav_control_all[i*24:(i+1)*24] = -1 * (vc[:, 30, 0] - vc[:, 32, 21] +
                                          vc[:, 31, 42] - vc[:, 37, 94] +
                                          vc[:, 39, 121] - vc[:, 34, 138] +
                                          vc[:, 35, 154] - vc[:, 30, 172])
        wav_all[i*24:(i+1)*24] = -1 * (vr[:, 30, 0] - vr[:, 32, 21] +
                                  vr[:, 31, 42] - vr[:, 37, 94] +
                                  vr[:, 39, 121] - vr[:, 34, 138] +
                                  vr[:, 35, 154] - vr[:, 30, 172])
        print('Done: ' + str(i+1) + ' out of ' + str(len(both)))
    circind = wav_all - wav_control_all
    for i in range(len(both)):
        circinds[0, i] = np.mean(np.concatenate((circind[24*i:24*i+2],
                                 circind[24*i+11:24*i+14],
                                 np.expand_dims(circind[24*i+23], axis=0)),
                                                axis=0), axis=0)
        circinds[1, i] = np.mean(np.concatenate((circind[24*i+5:24*i+8],
                                 circind[24*i+17:24*i+20]), axis=0), axis=0)
    anom = anom[:, ::-1, :]
    regmap_sig = np.zeros((2, 145, 192))
    regmap_sig[:] = np.NAN
    regmap_smoothed = np.zeros((2, 145, 192))
    regmap_smoothed_sig = np.zeros((2, 145, 192))
    regmap_smoothed_sig[:] = np.NAN
    for i in range(99):
        for j in range(192):
            for a in range(2):
                regmap_smoothed[a, i+23, j] = np.cov(anom[:, i+23, j], circinds[a, :])[0, 1]*3/(4* 6371**2 * 1.25 * 1.875 * (np.pi/180)**2)
                r = np.cov(anom[:, i+23, j], circinds[a, :])[0, 1]/np.sqrt(np.cov(anom[:, i+23, j], circinds[a, :])[1, 1]*np.cov(anom[:, i+23, j], circinds[a, :])[0, 0])
                t = r*np.sqrt((len(both)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i+23, j] = regmap_smoothed[a, i+23, j]
    return regmap_smoothed, regmap_smoothed_sig


def reconstruct(gto):
    from scipy import interpolate
    latR = np.linspace(90, -90, 73)
    lat = np.linspace(90, -90, 145)
    lonS = np.linspace(0, 358.125, 192)
    latS = np.linspace(88.541999816809453125, -88.541999816809453125, 94)
    sindices = np.zeros((36*3))
    windices = np.zeros((36*3))
    for i in range(36):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
    for i in range(36):
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    jetinds = np.zeros((2, 2, 36))
    u = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd.mon.mean.nc', 'uwnd')[:-25, 2, :]
    jetinds_monthly = np.zeros((2, 444))
    jetinds_monthly[0, :], jetinds_monthly[1, :] = jetindR(u, latR)
    wdata = jetinds_monthly[:, windices]
    sdata = jetinds_monthly[:, sindices]
    for y in range(36):
        jetinds[0, :, y] = np.mean(wdata[:, 3*y:3*(y+1)], axis=1)
        jetinds[1, :, y] = np.mean(sdata[:, 3*y:3*(y+1)], axis=1)
    jetind_d = np.zeros((2, 2, 36))
    for a in range(2):
        for b in range(2):
            jetind_d[a, b, :] = jetinds[a, b, :] - np.mean(jetinds[a, b, :])

    sst_n96 = np.zeros((444, 145, 192))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/skt.sfc.mon.mean.nc', 'skt')[:-12, :]
    for i in range(444):
        f = interpolate.interp2d(lonS, latS[::-1], sst[i, :])
        sst_n96[i, :] = f(lonS, lat)
    meshlat = np.zeros([np.ma.size(lat), 192])
    meshlat[:, :] = np.expand_dims(lat, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2

    jet_anom = np.zeros((2, 2, 36))
    wdataS = sst_n96[windices]
    sdataS = sst_n96[sindices]
    wdataSST = np.zeros((36, 145, 192))
    sdataSST = np.zeros((36, 145, 192))
    for y in range(36):
        wdataSST[y, :] = np.mean(wdataS[3*y:3*(y+1), :], axis=0)
        sdataSST[y, :] = np.mean(sdataS[3*y:3*(y+1), :], axis=0)
    # SST anom method 1
    #SSTclim = ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926.nc','SST_cpl')
    #wdataSST = wdataSST - np.mean(SSTclim[[0, 10, 11], ::-1, :], axis=0) - 273.15
    #sdataSST = sdataSST - np.mean(SSTclim[[5, 6, 7], ::-1, :], axis=0) - 273.15
    # SST anom method2
    wdataSST = wdataSST - np.mean(wdataSST, axis=0)
    sdataSST = sdataSST - np.mean(sdataSST, axis=0)

    wdataSST = wdataSST * meshlatweight
    sdataSST = sdataSST * meshlatweight
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    for b in range(2):
        for y in range(36):
            jet_anom[0, b, y] = np.nansum(wdataSST[y]*gto[0, b]*(1-lsm1))/20  # 2beta
            jet_anom[1, b, y] = np.nansum(sdataSST[y]*gto[1, b]*(1-lsm1))/20
    return jetind_d, jet_anom


def reconstruct_daily(gto, season='djf'):
    from scipy import interpolate
    import glob
    latR = np.linspace(90, -90, 73)
    lat = np.linspace(90, -90, 145)
    lonS = np.linspace(0, 358.125, 192)
    latS = np.linspace(88.541999816809453125, -88.541999816809453125, 94)
    if season == 'djf':
        months = 3
        sindices = np.zeros((36*3))
        windices = np.zeros((36*3))
        for i in range(36):
            sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    if season == 'ndj':
        months = 3
        sindices = np.zeros((36*3))
        windices = np.zeros((36*3))
        for i in range(36):
            sindices[3*i:3*(i+1)] = [16+12*i, 17+12*i, 18+12*i]
            windices[3*i:3*(i+1)] = [10+12*i, 11+12*i, 12+12*i]
    if season == 'son':
        months = 3
        sindices = np.zeros((36*3))
        windices = np.zeros((36*3))
        for i in range(36):
            sindices[3*i:3*(i+1)] = [14+12*i, 15+12*i, 16+12*i]
            windices[3*i:3*(i+1)] = [8+12*i, 9+12*i, 10+12*i]
    if season == 'nd':
        months = 2
        sindices = np.zeros((36*2))
        windices = np.zeros((36*2))
        for i in range(36):
            sindices[2*i:2*(i+1)] = [16+12*i, 17+12*i]
            windices[2*i:2*(i+1)] = [10+12*i, 11+12*i]
    if season == 'dj':
        months = 2
        sindices = np.zeros((36*2))
        windices = np.zeros((36*2))
        for i in range(36):
            sindices[2*i:2*(i+1)] = [17+12*i, 18+12*i]
            windices[2*i:2*(i+1)] = [11+12*i, 12+12*i]
    if season == 'n':
        months = 1
        sindices = np.zeros((36))
        windices = np.zeros((36))
        for i in range(36):
            sindices[i] = 16+12*i
            windices[i] = 10+12*i
    if season == 'd':
        months = 1
        sindices = np.zeros((36))
        windices = np.zeros((36))
        for i in range(36):
            sindices[i] = 17+12*i
            windices[i] = 11+12*i
    if season == 'on':
        months = 2
        sindices = np.zeros((36*2))
        windices = np.zeros((36*2))
        for i in range(36):
            sindices[2*i:2*(i+1)] = [15+12*i, 16+12*i]
            windices[2*i:2*(i+1)] = [9+12*i, 11+12*i]

    sindices = sindices.astype(int)
    windices = windices.astype(int)
    jetinds = np.zeros((2, 2, 36))
    jetinds_w = np.zeros((2, 3330))
    jetinds_s = np.zeros((2, 3404))
    a = glob.glob('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd_daily/years/*')
    for i in range(37):
        u_jf = ncread(a[i], 'uwnd')[:59, 2, :]
        u_jja = ncread(a[i], 'uwnd')[-214:-122, 2, :]
        u_d = ncread(a[i], 'uwnd')[-31:, 2, :]
        j_jf = jetindR(u_jf, latR)
        j_jja = jetindR(u_jja, latR)
        j_d = jetindR(u_d, latR)
        jetinds_w[:, 90*i:90*(i+1)] = np.concatenate((j_jf, j_d), axis=1)
        jetinds_s[:, 92*i:92*(i+1)] = j_jja
        #print(str(i))
    for y in range(36):
        jetinds[0, :, y] = np.mean(jetinds_w[:, 59+90*y:149+90*y], axis=1)
        jetinds[1, :, y] = np.mean(jetinds_s[:, 92+92*y:184+92*y], axis=1)
    jetind_d = np.zeros((2, 2, 36))
    for a in range(2):
        for b in range(2):
            jetind_d[a, b, :] = jetinds[a, b, :] - np.mean(jetinds[a, b, :])

    sst_n96 = np.zeros((444, 145, 192))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/skt.sfc.mon.mean.nc', 'skt')[:-12, :]
    for i in range(444):
        f = interpolate.interp2d(lonS, latS[::-1], sst[i, :])
        sst_n96[i, :] = f(lonS, lat)
    meshlat = np.zeros([np.ma.size(lat), 192])
    meshlat[:, :] = np.expand_dims(lat, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2

    jet_anom = np.zeros((2, 2, 36))
    wdataS = sst_n96[windices]
    sdataS = sst_n96[sindices]
    wdataSST = np.zeros((36, 145, 192))
    sdataSST = np.zeros((36, 145, 192))
    if months == 1:
        for y in range(36):
            wdataSST[y, :] = wdataS[months*y, :]
            sdataSST[y, :] = sdataS[months*y, :]
    else:
        for y in range(36):
            wdataSST[y, :] = np.mean(wdataS[months*y:months*(y+1), :], axis=0)
            sdataSST[y, :] = np.mean(sdataS[months*y:months*(y+1), :], axis=0)
    # SST anom method 1
    SSTclim = ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926.nc','SST_cpl')
    wdataSST = wdataSST - np.mean(SSTclim[[0, 10, 11], ::-1, :], axis=0) - 273.15
    sdataSST = sdataSST - np.mean(SSTclim[[5, 6, 7], ::-1, :], axis=0) - 273.15
    # SST anom method2
    #wdataSST = wdataSST - np.mean(wdataSST, axis=0)
    #sdataSST = sdataSST - np.mean(sdataSST, axis=0)
    wdataSST = wdataSST * meshlatweight
    sdataSST = sdataSST * meshlatweight
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    for b in range(2):
        for y in range(36):
            jet_anom[0, b, y] = np.nansum(wdataSST[y]*gto[0, b]*(1-lsm1))/20  # 2beta
            jet_anom[1, b, y] = np.nansum(sdataSST[y]*gto[1, b]*(1-lsm1))/20
    return jetind_d, jet_anom


def reconstruct_nao(gto, season='djf', rmask=np.ones((145, 192))):
    from scipy import interpolate
    lat = np.linspace(90, -90, 145)
    lonS = np.linspace(0, 358.125, 192)
    latS = np.linspace(88.541999816809453125, -88.541999816809453125, 94)
    sindices = np.zeros((36*3))
    windices = np.zeros((36*3))
    if season == 'djf':
        months = 3
        sindices = np.zeros((36*3))
        windices = np.zeros((36*3))
        for i in range(36):
            sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
            windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    if season == 'ndj':
        months = 3
        sindices = np.zeros((36*3))
        windices = np.zeros((36*3))
        for i in range(36):
            sindices[3*i:3*(i+1)] = [16+12*i, 17+12*i, 18+12*i]
            windices[3*i:3*(i+1)] = [10+12*i, 11+12*i, 12+12*i]
    if season == 'son':
        months = 3
        sindices = np.zeros((36*3))
        windices = np.zeros((36*3))
        for i in range(36):
            sindices[3*i:3*(i+1)] = [14+12*i, 15+12*i, 16+12*i]
            windices[3*i:3*(i+1)] = [8+12*i, 9+12*i, 10+12*i]
    if season == 'nd':
        months = 2
        sindices = np.zeros((36*2))
        windices = np.zeros((36*2))
        for i in range(36):
            sindices[2*i:2*(i+1)] = [16+12*i, 17+12*i]
            windices[2*i:2*(i+1)] = [10+12*i, 11+12*i]
    if season == 'dj':
        months = 2
        sindices = np.zeros((36*2))
        windices = np.zeros((36*2))
        for i in range(36):
            sindices[2*i:2*(i+1)] = [17+12*i, 18+12*i]
            windices[2*i:2*(i+1)] = [11+12*i, 12+12*i]
    if season == 'n':
        months = 1
        sindices = np.zeros((36))
        windices = np.zeros((36))
        for i in range(36):
            sindices[i] = 16+12*i
            windices[i] = 10+12*i
    if season == 'd':
        months = 1
        sindices = np.zeros((36))
        windices = np.zeros((36))
        for i in range(36):
            sindices[i] = 17+12*i
            windices[i] = 11+12*i
    if season == 'on':
        months = 2
        sindices = np.zeros((36*2))
        windices = np.zeros((36*2))
        for i in range(36):
            sindices[2*i:2*(i+1)] = [15+12*i, 16+12*i]
            windices[2*i:2*(i+1)] = [9+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    NAOinds = np.zeros((2, 36))
    slp_all = np.zeros((2, 12*37))
    slp = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/mslp.mon.mean.nc','mslp')[:-29]
    slp_all = slp[:, 21, 140], slp[:, 10, 135]
    slp_all = np.transpose(slp_all)
    for i in range(12):
        slp_all[i::12, :] = ((slp_all[i::12, :] -
                              np.mean(slp_all[i::12, :], axis=0)) /
                             np.std(slp_all[i::12, :], axis=0))
    NAOind = slp_all[:, 0] - slp_all[:, 1]
    for i in range(36):
        NAOinds[0, i] = np.mean(NAOind[12*i+11:12*i+14], axis=0)
        NAOinds[1, i] = np.mean(NAOind[12*i+17:12*i+20], axis=0)

    sst_n96 = np.zeros((444, 145, 192))
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/skt.sfc.mon.mean.nc', 'skt')[:-12, :]
    for i in range(444):
        f = interpolate.interp2d(lonS, latS[::-1], sst[i, :])
        sst_n96[i, :] = f(lonS, lat)
    meshlat = np.zeros([np.ma.size(lat), 192])
    meshlat[:, :] = np.expand_dims(lat, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1.25 * 1.875 * (np.pi/180)**2
    nao_anom = np.zeros((2, 36))   
    wdataS = sst_n96[windices]
    sdataS = sst_n96[sindices]
    wdataSST = np.zeros((36, 145, 192))
    sdataSST = np.zeros((36, 145, 192))
    if months == 1:
        for y in range(36):
            wdataSST[y, :] = wdataS[months*y, :]
            sdataSST[y, :] = sdataS[months*y, :]
    else:
        for y in range(36):
            wdataSST[y, :] = np.mean(wdataS[months*y:months*(y+1), :], axis=0)
            sdataSST[y, :] = np.mean(sdataS[months*y:months*(y+1), :], axis=0)
    # SST anom method 1
    #SSTclim = ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926.nc','SST_cpl')
    #wdataSST = wdataSST - np.mean(SSTclim[[0, 10, 11], ::-1, :], axis=0) - 273.15
    #sdataSST = sdataSST - np.mean(SSTclim[[5, 6, 7], ::-1, :], axis=0) - 273.15
    # SST anom method2
    #wdataSST = wdataSST - np.mean(wdataSST, axis=0)
    #sdataSST = sdataSST - np.mean(sdataSST, axis=0)

    wdataSST = wdataSST * meshlatweight
    sdataSST = sdataSST * meshlatweight
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
    return NAOinds, nao_anom


def scatter_sub(jet_act, jet_rec):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1980, 2016)
    for s in range(2):
        for j in range(2):
            jet_act_s = (jet_act[s, j] - np.mean(jet_act[s, j])) / np.std(jet_act[s, j])
            jet_rec_s = (jet_rec[s, j] - np.mean(jet_rec[s, j])) / np.std(jet_rec[s, j])
            axs[s, 1-j].plot(time, jet_act_s, color='k', label='Observed')
            axs[s, 1-j].plot(time, jet_rec_s, color='r', label='Reconstructed')
            axs[s, 1-j].set_ylim([-3, 3])
            axs[s, 1-j].set_xlim([1980, 2015])
            cc = np.corrcoef(jet_act_s, jet_rec_s)[0, 1]
            axs[s, 1-j].annotate('$\mathregular{r}=$%.2f' % cc, xy=(2010.5, 2.6), fontsize=16)
            axs[s, 1-j].set_xlabel('Year', fontsize=16)
            axs[s, 1-j].legend(loc=2, fontsize=16)
            if j == 0:
                axs[s, 1-j].set_ylabel('Jet speed anomaly (ms$\mathregular{^{-1}}$)', fontsize=16)
                if s == 0:
                    axs[s, 1-j].set_title('Speed', fontsize=20)
            if j == 1:
                axs[s, 1-j].set_ylabel('Jet position anomaly (deg)', fontsize=16)
                if s == 0:
                    axs[s, 1-j].set_title('Latitude', fontsize=16)
            axs[s, 1-j].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def scatter_subnao(nao_act, nao_rec):
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 1, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    time = np.arange(1980, 2016)
    for s in range(2):
        jet_act_s = (nao_act[s] - np.mean(nao_act[s])) / np.std(nao_act[s])
        jet_rec_s = (nao_rec[s] - np.mean(nao_rec[s])) / np.std(nao_rec[s])
        axs[s].plot(time, jet_act_s, color='k', label='Observed')
        axs[s].plot(time, jet_rec_s, color='r', label='Reconstructed')
        axs[s].set_ylim([-3, 3])
        axs[s].set_xlim([1980, 2015])
        cc = np.corrcoef(jet_act_s, jet_rec_s)[0, 1]
        axs[s].annotate('$\mathregular{r}=$%.2f' % cc, xy=(2012.5, 2.6),
                        fontsize=16)
        axs[s].set_xlabel('Year', fontsize=16)
        axs[s].legend(loc=2, fontsize=16)
        if s == 0:
            axs[s].set_title('NAO', fontsize=16)
        axs[s].tick_params(axis='both', which='major', labelsize=16)
    rows = ['Winter', 'Summer']
    for ax, row in zip(axs, rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-55, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=20)
    for ax, lab in zip(axs, ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=20, fontweight='bold')
    plt.subplots_adjust(hspace=0.25, wspace=0.15, top=.95, bottom=0.075, left=.075,
                        right=.95)


def mapplot(plotdata, lat, lon, mx=2, mask='no', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap, shiftgrid
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        plotdata = np.ma.masked_array(plotdata, lsm)
    long = np.concatenate((lon[95:]-360, lon[:97]))
    plotdata = np.concatenate((plotdata[:, 95:], plotdata[:, :97]), axis=1)
    plotdata1, lon = shiftgrid(180., plotdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)
    plt.figure()
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[239:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
    for i in range(np.ma.size(lat)):
        for j in range(np.ma.size(lon)):
            if plotdata[i, j] == 'masked':
                plotdata[i, j] = np.nan
    plot = m.contourf(x, y, plotdata, ctrs,
                      cmap=newcmap, vmin=caxismin, vmax=caxismax,
                      extend='both')
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def mapplot60(plotdata, plotdata2, lat, lon, mask='no', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        for b in range(2):
            if b == 0:
                mx = 1
            else:
                mx = 1
            plotdata1 = np.concatenate((plotdata[a, 1-b, :],
                                        np.expand_dims(plotdata[a, 1-b, :, 0],
                                        axis=1)), axis=1)
            plotdata3 = np.concatenate((plotdata2[a, 1-b, :],
                                        np.expand_dims(plotdata2[a, 1-b, :, 0],
                                        axis=1)), axis=1)
            if mask == 'yes':
                lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
                lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                     axis=1)
                plotdata1 = np.ma.masked_array(plotdata1, lsm)
                plotdata3 = np.ma.masked_array(plotdata3, lsm)
            meshlon, meshlat = np.meshgrid(long, lat)
            ax1 = axs[a, b]
            m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                        llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
            m.drawcoastlines()
            m.drawmapboundary(linewidth=2)
            x, y = m(meshlon, meshlat)
            mycmap2 = plt.cm.YlOrRd(np.arange(256))
            mycmap1 = plt.cm.Blues_r(np.arange(256))
            my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
            #my_cmap[239:274, :] = 1
            newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                                   my_cmap)
            caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
            plot = m.contourf(x, y, plotdata1, ctrs,
                              cmap=newcmap, vmin=caxismin, vmax=caxismax,
                              extend='both')
            m.contour(x, y, plotdata3, ctrs, colors='k')
            ax1.set_ylim(-60, 60)
            m.drawparallels(np.arange(-60, 90, 30),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0, 390, 30),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-60., 90., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if a == 0 and b == 0:
                ax1.set_ylabel('Winter', fontsize=16, labelpad=25)
                ax1.set_title('Latitude', fontsize=16, y=1.08)
            if a == 1 and b == 0:
                c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
                ax1.set_ylabel('Summer', fontsize=16, labelpad=25)
                c.set_label(label='Poleward jet latitude shift ($^\circ$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
            if a == 0 and b == 1:
                ax1.set_title('Speed', fontsize=16, y=1.08)
            if a == 1 and b == 1:
                c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
                c.set_label(label='Jet speed increase (ms$\mathregular{^{-1}}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    plt.subplots_adjust(hspace=0, wspace=0.1, top=.95, bottom=0.05, left=.05,
                        right=.95)

    plt.show()


def mapplot60_nao(plotdata_sig, plotdata, lat, lon, mx=2, mask='no', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 1, facecolor='w',
                            edgecolor='k', linewidth=2)
    ctrs2 = np.concatenate((np.arange(-2, 0, .25), np.arange(.25, 2.25, .25)))
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        ax1 = axs[a]
        plotdata1 = np.concatenate((plotdata_sig[a],
                                    np.expand_dims(plotdata[a, :, 0],
                                    axis=1)), axis=1)
        plotdata2 = np.concatenate((plotdata[a],
                                    np.expand_dims(plotdata[a, :, 0],
                                    axis=1)), axis=1)
        if mask == 'yes':
            lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
            lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                 axis=1)
            plotdata1 = np.ma.masked_array(plotdata1, lsm)
            plotdata2 = np.ma.masked_array(plotdata2, lsm)
        meshlon, meshlat = np.meshgrid(long, lat)
        m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                    llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
        m.drawcoastlines()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        mycmap2 = plt.cm.YlOrRd(np.arange(256))
        mycmap1 = plt.cm.Blues_r(np.arange(256))
        my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
        #my_cmap[238:275, :] = 1
        newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                               my_cmap)
        caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
        plot = m.contourf(x, y, plotdata1, ctrs,
                          cmap=newcmap, vmin=caxismin, vmax=caxismax,
                          extend='both')
        m.contour(x, y, plotdata2, ctrs, colors='k')
        ax1.set_ylim(-60, 60)
        m.drawparallels(np.arange(-60., 90., 30.),
                        labels=[True, False, False, True], linewidth=0)
        m.drawmeridians(np.arange(0., 390., 30.),
                        labels=[True, False, False, True], linewidth=0)
        ax1.set_yticks(np.arange(-60., 90., 30.))
        ax1.set_xticks(np.arange(0., 390., 30.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
        ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)
        if a == 0:
            ax1.set_ylabel('Winter', fontsize=16, labelpad=20)
        else:
            ax1.set_ylabel('Summer', fontsize=16, labelpad=20)
    c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                     spacing='proportional', aspect=50)
    c.set_label(label='NAO sensitivity ([$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    #c.set_label(label='Circulation index sensitivity (ms$^{-1}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    plt.subplots_adjust(hspace=0, wspace=0.1, top=.95, bottom=0.05, left=.05,
                        right=.95)
    plt.show()













