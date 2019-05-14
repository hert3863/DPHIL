# -*- coding: utf-8 -*-
"""
Created on Thu Jun  2 13:20:23 2016

@author: bakerh

NAME
    GCM Analysis
PURPOSE
    Script to analyse results of iGCM run
"""
import numpy as np
from netCDF4 import Dataset
import os
from scipy import stats

'''
MUST PERFORM DAILY ANALYSIS OUTPUT BEFORE RUNNING
'''


def I2Kfull(data):
    '''
    Inserts 2Kruna levels in to 2Krun levels
    '''
    data1 = np.copy(data)
    data1[35:69] = data[171:205]
    data1[69:103] = data[35:69]
    data1[103:137] = data[205:239]
    data1[137:171] = data[69:103]
    data1[171:205] = data[239:273]
    data1[205:239] = data[103:137]
    data1[239:273] = data[273:307]
    data1[273:307] = data[137:171]
    return data1


def indiv(experiment, run, invariable):
    '''
    outputs numpy array combining netcdfs for a run in slices

    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    import numpy as np
    import os
    from netCDF4 import Dataset
    yrs = []
    for d in os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        experiment + '/run' + str(run) + '/history/'):
        if d.find('daily.nc') == -1 and d.find('.', 0, 1) == -1:
            yrs.append(d)
    variable = np.zeros((37, 64))
    variable = np.zeros((len(yrs), 37, 64))
    for i in range(len(yrs)):
        nc_f = ('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                experiment + '/run' + str(run) + '/history/' + yrs[i])
        nc_fid = Dataset(nc_f, 'r')
        data = nc_fid.variables[invariable][:]
        variable[i, :, :] = data
    return variable


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


def combine(experiment, run, invariable):
    '''
    outputs numpy array combining netcdfs for a run

    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    variable = np.zeros((37, 64))
    months = []
    for d in os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        experiment + '/' + run + '/history/'):
            if d.find('daily.nc') == -1:
                months.append(d)
    reqmonths = months[12:]  # first 12 months spinup
    for i in reqmonths:
        nc_f = ('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                experiment + '/' + run + '/history/' + i)
        nc_fid = Dataset(nc_f, 'r')
        data = nc_fid.variables[invariable][:]
        variable = data + variable
    variable = variable / 72  # 72 months of data
    return variable


def multiread(experiment, invariable):
    '''
    outputs numpy array of all runs in experiment

    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    s = len(os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                       experiment))
    variables = np.zeros((s, 37, 64))
    for i in range(len(os.listdir('/network/aopp/hera/mad/' +
                                  'bakerh/data/FMS/output/' + experiment))):
        variables[i, :, :] = combine(experiment, 'run' + str(i+1), invariable)
    return variables


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
    u850 = u[:, 31, :]
    indices = np.zeros((np.ma.size(u850, axis=0), 2, 2))

    for t in range(np.ma.size(u850, axis=0)):
        indices[t, 0, 0] = np.amax(u850[t, 0:32])
        indices[t, 0, 1] = lat[np.argmax(u850[t, 0:32])]
        indices[t, 1, 0] = np.amax(u850[t, 39:])
        indices[t, 1, 1] = lat[np.argmax(u850[t, 39:])+39]
    return indices


def jetindicesdaily(u850, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    iGCM output

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
    speeds = np.zeros((np.ma.size(u850, axis=0), 2))
    lats = np.zeros((np.ma.size(u850, axis=0), 2))
    for t in range(np.ma.size(u850, axis=0)):
        wd, ww, sd, sw = jetdist(u850[t, :], lat)
        wm, sm = np.argmax(ww), np.argmax(sw)
        wp = np.polyfit(wd[wm-1:wm+2], ww[wm-1:wm+2], 2)
        sp = np.polyfit(sd[sm-1:sm+2], sw[sm-1:sm+2], 2)
        '''
        x = Symbol('x', real=True)
        f = wp[0] * x**2 + wp[1] * x + wp[2]
        fprime = f.diff(x)
        xmax = solve(fprime, x)
        ymax = wp[0] * xmax[0]**2 + wp[1] * xmax[0] + wp[2]
        a = Symbol('a', real=True)
        b = sp[0] * a**2 + sp[1] * a + sp[2]
        bprime = b.diff(a)
        amax = solve(bprime, a)
        bmax = sp[0] * amax[0]**2 + sp[1] * amax[0] + sp[2]
        '''
        def f(x): return wp[0] * x**2 + wp[1] * x + wp[2]
        xmax = optimize.fmin(lambda x: -f(x), wm, disp=False)
        ymax = f(xmax)

        def g(a): return sp[0] * a**2 + sp[1] * a + sp[2]
        amax = optimize.fmin(lambda a: -g(a), sm, disp=False)
        bmax = g(amax)
        speeds[t, 0] = ymax
        lats[t, 0] = xmax[0]
        speeds[t, 1] = bmax
        lats[t, 1] = amax[0]
    return speeds, lats


def cplot(plotdata, sigma, lat, title=''):
    """
    Plots input grid of quantity at sigma and lat coords

    Parameters
    ----------
    plotdata: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    plt.figure()
    my_cmap = plt.cm.RdYlBu_r(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    caxismin, caxismax, ctrs = colourscale(plotdata)
    plot = plt.contourf(meshlat, meshsigma, plotdata, ctrs,
                        cmap=newcmap, vmin=caxismin, vmax=caxismax)
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    plt.gca().invert_yaxis()
    plt.yscale('log')
    plt.xlim([-90, 90])
    plt.xticks(np.arange(-90, 105, 15))
    plt.title(title, y=1.08, fontsize=30)
    plt.show()
    # plt.clabel(plot, inline=0, fontsize=10, fmt='%1.1f', inline_spacing=1)
    plt.show()


def colourscale(plotdata):
    """
    Takes data being plotted and normalises the colourscale between largest
    data value and it's negative multiple

    Parameters
    ----------
    plotdata: array
        data being plotted

    Returns
    -------
    caxismax: int
        max magnitude data value
    caxismin: int
        negative of max mag data value
    ctrs: array
        gradated colour scale contour array
    """
    M = np.nanmax(plotdata)
    m = np.nanmin(plotdata)
    if M >= abs(m):
        ctrs1 = np.arange(-M, 0, .1*M)
        ctrs2 = np.arange(0.1*M, 1.09*M, .1*M)
        ctrs = np.concatenate((ctrs1, ctrs2))
        caxismin = -M
        caxismax = M
    else:
        m = -m
        ctrs1 = np.arange(-m, 0, .1*m)
        ctrs2 = np.arange(0.1*m, 1.09*m, .1*m)
        ctrs = np.concatenate((ctrs1, ctrs2))
        caxismin = -m
        caxismax = m
    # function will not work if there exist no positive max or negative min
    return caxismin, caxismax, ctrs


def jetdist(u850, lat):
    """
    Determines the winter and summer jet distributions

    Parameters
    ----------
    u850: array
        wind data at u850
    lat: array
        latitudes of data

    Returns
    -------
    wdist: array
        winter lat distribution
    wweights: array
        speeds at each lat
    sdist: array
        summer lat distribution
    sweights: array
        speeds at each lat
    """
    wwmx = np.max(u850[0:32])
    windex = np.argmax(u850[0:32])
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

    windex = np.argmax(u850[0:32])
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

    swmx = np.max(u850[32:])
    sindex = np.argmax(u850[32:]) + 32
    sdmax = lat[sindex]
    sdist = [sdmax]
    sweights = [swmx]

    sindex -= 1
    while u850[sindex+1] >= u850[sindex] and u850[sindex+1] >= 0:
        if u850[sindex] < 0:
            lat_0 = lat[sindex] - u850[sindex] * ((lat[sindex+1]-lat[sindex]) /
                                                  (u850[sindex+1] -
                                                   u850[sindex]))
            sdist.insert(0, lat_0)
            sweights.insert(0, 0)
            sindex -= 1
        else:
            sdist.insert(0, lat[sindex])
            sweights.insert(0, u850[sindex])
            sindex -= 1

    sindex = np.argmax(u850[32:]) + 32
    sindex += 1
    while sindex < 64 and u850[sindex-1] >= u850[sindex] and u850[sindex-1] >= 0:
        if u850[sindex] < 0:
            lat_0 = lat[sindex-1] - u850[sindex-1] * ((lat[sindex] -
                                                      lat[sindex-1]) /
                                                      (u850[sindex] -
                                                       u850[sindex-1]))
            sdist.append(lat_0)
            sweights.append(0)
            sindex += 1
        else:
            sdist.append(lat[sindex])
            sweights.append(u850[sindex])
            sindex += 1

    return wdist, wweights, sdist, sweights


def recast(data):
    data = np.reshape(data, (2, 2, 9, 34))
    return data


def sigtest(series1, series2):
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
    var1 = np.var(series1)
    var2 = np.var(series2)
    rho1 = np.corrcoef(series1[:-1], series1[1:])
    rho1 = rho1[1, 0]
    rho2 = np.corrcoef(series2[:-1], series2[1:])
    rho2 = rho2[1, 0]
    neff1 = np.ma.size(series1) * (1-rho1) / (1+rho1)
    neff2 = np.ma.size(series2) * (1-rho2) / (1+rho2)
    z = ((np.mean(series1) - np.mean(series2)) /
         (np.sqrt((var1 / neff1) + (var2 / neff2))))
    p = 1 - st.norm.cdf(np.abs(z))
    return p, z


def sigtest_monthly(series1, series2):
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
    var1 = np.var(series1)
    var2 = np.var(series2)
    z = ((np.mean(series1) - np.mean(series2)) /
         (np.sqrt((var1 / len(series1)) + (var2 / len(series2)))))
    p = 1 - st.norm.cdf(np.abs(z))
    return p, z


def maincoarse():
    sigma = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                   'coarserun27.5.16/run0/history/day3630h00.nc', 'sigma')
    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 'coarserun27.5.16/run0/history/day3630h00.nc', 'lat')
    '''
    speeds, lats = jetindicesmean(u, lat)
    d_speeds = np.zeros((70, 2))
    d_lats = np.zeros((70, 2))
    for a in range(70):
        for b in range(2):
            d_speeds[a, b] = speeds[a+1, b] - speeds[0, b]
            d_lats[a, b] = lats[a+1, b] - lats[0, b]
    d_lats = recast(d_lats)
    d_speeds = recast(d_speeds)
    # change sign of negative latitudes
    d_lats[:, :, :6] = d_lats[:, :, :6] * -1
    '''
    gridlat = np.zeros((14))
    gridlat[0], gridlat[1], gridlat[2] = lat[1], lat[6], lat[11]
    gridlat[3], gridlat[4], gridlat[5] = lat[16], lat[21], lat[26]
    gridlat[6], gridlat[7], gridlat[8] = lat[31], lat[32], lat[37]
    gridlat[9], gridlat[10], gridlat[11] = lat[42], lat[47], lat[52]
    gridlat[12], gridlat[13] = lat[57], lat[62]
    gridsigma = np.zeros((5))
    gridsigma[0], gridsigma[1], gridsigma[2] = sigma[12], sigma[18], sigma[24]
    gridsigma[3], gridsigma[4] = sigma[30], sigma[36]

    daily = np.zeros((71, 2160, 2, 2))
    u850 = np.zeros((71, 2160, 64))
    # run, day, winter/summer, speed/lat
    for i in range(3):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'coarserun27.5.16/run' + str(i) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)

    for i in range(4, 8):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'coarserun27.5.16/run' + str(i) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)

    for i in range(9, 13):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'coarserun27.5.16/run' + str(i) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)

    for i in range(14, 71):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'coarserun27.5.16/run' + str(i) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)

    means = np.mean(daily, axis=1)
    sd = np.std(daily, axis=1)
    skewn = stats.skew(daily, axis=1)
    kurt = stats.kurtosis(daily, axis=1, fisher=True)
    d_means = np.zeros((2, 2, 70))
    d_sd = np.zeros((2, 2, 70))
    d_skewn = np.zeros((2, 2, 70))
    d_kurt = np.zeros((2, 2, 70))
    for a in range(70):
        for b in range(2):
            for c in range(2):
                d_means[b, c, a] = means[a+1, b, c] - means[0, b, c]
                d_sd[b, c, a] = sd[a+1, b, c] - sd[0, b, c]
                d_skewn[b, c, a] = skewn[a+1, b, c] - skewn[0, b, c]
                d_kurt[b, c, a] = kurt[a+1, b, c] - kurt[0, b, c]
    d_means = recast(d_means)
    d_sd = recast(d_sd)
    d_skewn = recast(d_skewn)
    d_kurt = recast(d_kurt)
    # change sign of negative latitudes
    d_means[0, 1, :, :] = d_means[0, 1, :, :] * -1
    d_skewn[0, 1, :, :] = d_skewn[0, 1, :, :] * -1

    jetmeans = np.zeros((71, 2))
    jetsdmeans = np.zeros((71, 2))
    jetskewmeans = np.zeros((71, 2))
    jetkurtmeans = np.zeros((71, 2))
    for i in range(71):
        jetmean = np.zeros((2160, 2))
        jetvari = np.zeros((2160, 2))
        jetskewn = np.zeros((2160, 2))
        jetkurt = np.zeros((2160, 2))
        for j in range(2160):
            wdist, wweights, sdist, sweights = jetdist(u850[i, j, :], lat)
            wdis = stats.rv_discrete(values=(wdist, wweights/np.sum(wweights)))
            sdis = stats.rv_discrete(values=(sdist, sweights/np.sum(sweights)))
            jetmean[j, 0], jetvari[j, 0], jetskewn[j, 0], jetkurt[j, 0] = wdis.stats(moments='mvsk')
            jetmean[j, 1], jetvari[j, 1], jetskewn[j, 1], jetkurt[j, 1] = sdis.stats(moments='mvsk')
        jetmeans[i, :] = np.mean(jetmean, axis=0)
        jetsdmeans[i, :] = np.mean(np.sqrt(jetvari), axis=0)
        jetskewmeans[i, :] = np.mean(jetskewn, axis=0)
        jetkurtmeans[i, :] = np.mean(jetkurt, axis=0)

    d_jetmean = np.zeros((2, 70))
    d_jetsdmean = np.zeros((2, 70))
    d_jetskewmean = np.zeros((2, 70))
    d_jetkurtmean = np.zeros((2, 70))
    for a in range(70):
        for b in range(2):
                d_jetmean[b, a] = jetmeans[a+1, b] - jetmeans[0, b]
                d_jetsdmean[b, a] = jetsdmeans[a+1, b] - jetsdmeans[0, b]
                d_jetskewmean[b, a] = jetskewmeans[a+1, b] - jetskewmeans[0, b]
                d_jetkurtmean[b, a] = jetkurtmeans[a+1, b] - jetkurtmeans[0, b]
    d_jetmean = np.reshape(d_jetmean, (2, 5, 14))
    d_jetsdmean = np.reshape(d_jetsdmean, (2, 5, 14))
    d_jetskewmean = np.reshape(d_jetskewmean, (2, 5, 14))
    d_jetkurtmean = np.reshape(d_jetkurtmean, (2, 5, 14))
    d_jetmean[0, :] = d_jetmean[0, :] * -1
    d_jetskewmean[0, :] = d_jetskewmean[0, :] * -1  # change sign of winter
    return d_means, d_sd, d_skewn, d_kurt, d_jetmean, d_jetsdmean, d_jetskewmean, d_jetkurtmean


def main2K():
    sigma = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                   '2Krun/run0/history/day3630h00.nc', 'sigma')
    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')
    '''
    speeds, lats = jetindicesmean(u, lat)
    d_speeds = np.zeros((70, 2))
    d_lats = np.zeros((70, 2))
    for a in range(70):
        for b in range(2):
            d_speeds[a, b] = speeds[a+1, b] - speeds[0, b]
            d_lats[a, b] = lats[a+1, b] - lats[0, b]
    d_lats = recast(d_lats)
    d_speeds = recast(d_speeds)
    # change sign of negative latitudes
    d_lats[:, :, :6] = d_lats[:, :, :6] * -1
    '''
    gridlat = np.concatenate((lat[0:32:2], lat[31:33], lat[33::2]))
    gridsigma = np.zeros((9))
    gridsigma[0], gridsigma[1], gridsigma[2] = sigma[36], sigma[30], sigma[26]
    gridsigma[3], gridsigma[4] = sigma[22], sigma[16]
    gridsigma[5], gridsigma[6], gridsigma[7] = sigma[33], sigma[28], sigma[24]
    gridsigma[8] = sigma[19]

    daily = np.zeros((307, 2160, 2, 2))
    # u850 = np.zeros((307, 2160, 64))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        '2Krun/run' + str(i) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        # u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)
        print('Completed polyfit for run' + str(i))

    means = np.mean(daily, axis=1)
    sd = np.std(daily, axis=1)
    skewn = stats.skew(daily, axis=1)
    kurt = stats.kurtosis(daily, axis=1, fisher=True)
    d_means = np.zeros((2, 2, 306))
    d_sd = np.zeros((2, 2, 306))
    d_skewn = np.zeros((2, 2, 306))
    d_kurt = np.zeros((2, 2, 306))
    for a in range(306):
        for b in range(2):
            for c in range(2):
                d_means[b, c, a] = means[a+1, b, c] - means[0, b, c]
                d_sd[b, c, a] = sd[a+1, b, c] - sd[0, b, c]
                d_skewn[b, c, a] = skewn[a+1, b, c] - skewn[0, b, c]
                d_kurt[b, c, a] = kurt[a+1, b, c] - kurt[0, b, c]
    d_means = recast(d_means)
    d_sd = recast(d_sd)
    d_skewn = recast(d_skewn)
    d_kurt = recast(d_kurt)
    # change sign of negative latitudes
    d_means[0, 1, :, :] = d_means[0, 1, :, :] * -1
    d_skewn[0, 1, :, :] = d_skewn[0, 1, :, :] * -1
    '''
    jetmeans = np.zeros((171, 2))
    jetsdmeans = np.zeros((171, 2))
    jetskewmeans = np.zeros((171, 2))
    jetkurtmeans = np.zeros((171, 2))
    for i in range(171):
        jetmean = np.zeros((2160, 2))
        jetvari = np.zeros((2160, 2))
        jetskewn = np.zeros((2160, 2))
        jetkurt = np.zeros((2160, 2))
        for j in range(2160):
            wdist, wweights, sdist, sweights = jetdist(u850[i, j, :], lat)
            wdis = stats.rv_discrete(values=(wdist, wweights/np.sum(wweights)))
            sdis = stats.rv_discrete(values=(sdist, sweights/np.sum(sweights)))
            jetmean[j, 0], jetvari[j, 0], jetskewn[j, 0], jetkurt[j, 0] = wdis.stats(moments='mvsk')
            jetmean[j, 1], jetvari[j, 1], jetskewn[j, 1], jetkurt[j, 1] = sdis.stats(moments='mvsk')
        jetmeans[i, :] = np.mean(jetmean, axis=0)
        jetsdmeans[i, :] = np.mean(np.sqrt(jetvari), axis=0)
        jetskewmeans[i, :] = np.mean(jetskewn, axis=0)
        jetkurtmeans[i, :] = np.mean(jetkurt, axis=0)
        print('Completed jetdist for run' + str(i))

    d_jetmean = np.zeros((2, 170))
    d_jetsdmean = np.zeros((2, 170))
    d_jetskewmean = np.zeros((2, 170))
    d_jetkurtmean = np.zeros((2, 170))
    for a in range(170):
        for b in range(2):
                d_jetmean[b, a] = jetmeans[a+1, b] - jetmeans[0, b]
                d_jetsdmean[b, a] = jetsdmeans[a+1, b] - jetsdmeans[0, b]
                d_jetskewmean[b, a] = jetskewmeans[a+1, b] - jetskewmeans[0, b]
                d_jetkurtmean[b, a] = jetkurtmeans[a+1, b] - jetkurtmeans[0, b]
    d_jetmean = np.reshape(d_jetmean, (2, 5, 34))
    d_jetsdmean = np.reshape(d_jetsdmean, (2, 5, 34))
    d_jetskewmean = np.reshape(d_jetskewmean, (2, 5, 34))
    d_jetkurtmean = np.reshape(d_jetkurtmean, (2, 5, 34))
    d_jetmean[0, :] = d_jetmean[0, :] * -1
    d_jetskewmean[0, :] = d_jetskewmean[0, :] * -1  # change sign of winter
    '''
    return daily, d_means, d_sd, d_skewn, d_kurt  # , d_jetmean, d_jetsdmean, d_jetskewmean, d_jetkurtmean


def main2Kst(variable):
    sigma = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                   '2Krun/run0/history/day3630h00.nc', 'sigma')
    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')

    gridlat = np.concatenate((lat[0:32:2], lat[31:33], lat[33::2]))
    gridsigma = np.zeros((9))
    gridsigma[0], gridsigma[1], gridsigma[2] = sigma[36], sigma[30], sigma[26]
    gridsigma[3], gridsigma[4] = sigma[22], sigma[16]
    gridsigma[5], gridsigma[6], gridsigma[7] = sigma[33], sigma[28], sigma[24]
    gridsigma[8] = sigma[19]

    daily = np.zeros((307, 72, 2, 2))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        udaily = indiv('2Krun', str(i), variable)[12:, :, :]
        #  for t in range(72):
        #  udaily[t,:] = div(udaily[t,:], sigma, lat)
        u850 = udaily[:, 31, :]  # 18 for uv
        mag = np.zeros((np.ma.size(u850, axis=0), 2))
        lats = np.zeros((np.ma.size(u850, axis=0), 2))
        for t in range(72):
            mag[t, 0] = np.amin(u850[t, 0:32])
            lats[t, 0] = lat[np.argmin(u850[t, 0:32])]
            mag[t, 1] = np.amax(u850[t, 39:])  # min for uv
            lats[t, 1] = lat[np.argmax(u850[t, 39:])+39]  # min for uv
        daily[i, :, :, 0], daily[i, :, :, 1] = mag, lats
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    sd = np.std(daily, axis=1)
    skewn = stats.skew(daily, axis=1)
    kurt = stats.kurtosis(daily, axis=1, fisher=True)
    d_means = np.zeros((2, 2, 306))
    d_sd = np.zeros((2, 2, 306))
    d_skewn = np.zeros((2, 2, 306))
    d_kurt = np.zeros((2, 2, 306))
    for a in range(306):
        for b in range(2):
            for c in range(2):
                d_means[b, c, a] = means[a+1, b, c] - means[0, b, c]
                d_sd[b, c, a] = sd[a+1, b, c] - sd[0, b, c]
                d_skewn[b, c, a] = skewn[a+1, b, c] - skewn[0, b, c]
                d_kurt[b, c, a] = kurt[a+1, b, c] - kurt[0, b, c]
    d_means = recast(d_means)
    d_sd = recast(d_sd)
    d_skewn = recast(d_skewn)
    d_kurt = recast(d_kurt)
    # change sign of negative latitudes
    d_means[0, 1, :, :] = d_means[0, 1, :, :] * -1
    d_skewn[0, 1, :, :] = d_skewn[0, 1, :, :] * -1
    # change of sign for negative mag (need to do summer speeds too for uv)
    d_means[0, 0, :, :] = d_means[0, 0, :, :] * -1

    return daily, d_means, #  d_sd, d_skewn, d_kurt, d_jetmean, d_jetsdmean, d_jetskewmean, d_jetkurtmean


def paper1_q0(jetindcontrol, lat):
    daily = np.zeros((20, 2160, 2, 2))
    # u850 = np.zeros((307, 2160, 64))
    # run, day, winter/summer, speed/lat
    for i in range(20):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'paper1_q0_1/run' + str(i+1) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        # u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)
        print('Completed polyfit for run' + str(i+1))

    means = np.mean(daily, axis=1)
    sd = np.std(daily, axis=1)
    d_means = np.zeros((20, 2, 2))
    for a in range(20):
        for b in range(2):
            for c in range(2):
                d_means[a, b, c] = means[a, b, c] - jetindcontrol[b, c]

    # change sign of negative latitudes
    d_means[:, 0, 1] = d_means[:, 0, 1] * -1

    return daily, d_means, sd


def paper1_q0b(jetindcontrol, lat):
    daily = np.zeros((4, 2160, 2, 2))
    # u850 = np.zeros((307, 2160, 64))
    # run, day, winter/summer, speed/lat
    for i in range(4):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'paper1_q0_1b/run' + str(17+i) + '/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        # u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)
        print('Completed polyfit for run' + str(i+1))

    means = np.mean(daily, axis=1)
    sd = np.std(daily, axis=1)
    d_means = np.zeros((4, 2, 2))
    for a in range(4):
        for b in range(2):
            for c in range(2):
                d_means[a, b, c] = means[a, b, c] - jetindcontrol[b, c]

    # change sign of negative latitudes
    d_means[:, 0, 1] = d_means[:, 0, 1] * -1

    return daily, d_means, sd


def paper1_qadd(jetindcontrol, lat):
    daily = np.zeros((11, 2160, 2, 2))
    # u850 = np.zeros((307, 2160, 64))
    # run, day, winter/summer, speed/lat
    for i in range(11):
        udaily = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        'paper1_qadd/run' + str(i+1) + '/history/daily.nc',
                        'uwnd')
        udaily1 = udaily[360:, 2, :]
        # u850[i, :, :] = udaily1
        daily[i, :, :, 0], daily[i, :, :, 1] = jetindicesdaily(udaily1, lat)
        print('Completed polyfit for run' + str(i))

    means = np.mean(daily, axis=1)
    sd = np.std(daily, axis=1)
    d_means = np.zeros((11, 2, 2))
    for a in range(11):
        for b in range(2):
            for c in range(2):
                d_means[a, b, c] = means[a, b, c] - jetindcontrol[b, c]

    # change sign of negative latitudes
    d_means[:, 0, 1] = d_means[:, 0, 1] * -1

    return daily, d_means, sd


def main2KHad():
    def getHadEx(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(50.0, lat)
            strfcn1 = strfcn[:, idx_min:idx_max+1] * -1.0
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            strfcn1 = strfcn[:, idx_min:idx_max+1]
            lat1 = lat[idx_min:idx_max+1]

        sfmin = np.min(strfcn1, axis=0)  # Find minimum val strfcn in column
        imin = np.argmin(strfcn1, axis=0)  # Find the sig index of the minimum in column
        a = np.argmin(sfmin)            # Find the lat index of global minimum
        lev_sfmin = imin[a]            # Find the sig index of global minimum
        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[lev_sfmin, a]) == np.sign(strfcn1[lev_sfmin, aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[lev_sfmin, aa+1],
                                               strfcn1[lev_sfmin, aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < sfmin.size and np.sign(strfcn1[lev_sfmin, a]) == np.sign(strfcn1[lev_sfmin, aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[lev_sfmin, aa-1],
                                               strfcn1[lev_sfmin, aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))

    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')

    daily = np.zeros((307, 72, 2))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        hadext = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadEx(strfcn[t, :], lat, hem='NH')
            hadext[t, 1] = getHadEx(strfcn[t, :], lat, hem='SH')

        daily[i, :, 1], daily[i, :, 0] = hadext[:, 0], hadext[:, 1]
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    means = I2Kfull(means)
    d_means = np.zeros((2, 306))
    #  d_sd = np.zeros((2, 306))
    for a in range(306):
        for b in range(2):
            d_means[b, a] = means[a+1, b] - means[0, b]

    d_means = np.reshape(d_means, (2, 9, 34))
    # change sign of negative latitudes
    d_means[0, :] = d_means[0, :] * -1

    return daily, means, d_means


def main2KHad_vertav():
    def getHadEx(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(50.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
            lat1 = lat[idx_min:idx_max+1]

        a = np.argmin(strfcn1)           # Find the lat index of global minimum

        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))

    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')

    daily = np.zeros((307, 72, 2))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        hadext = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadEx(strfcn[t, :], lat, hem='NH')
            hadext[t, 1] = getHadEx(strfcn[t, :], lat, hem='SH')

        daily[i, :, 1], daily[i, :, 0] = hadext[:, 0], hadext[:, 1]
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    stand = np.std(daily, axis=1)
    means = I2Kfull(means)
    stand = I2Kfull(stand)
    d_means = np.zeros((2, 306))
    d_std = np.zeros((2, 306))
    for a in range(306):
        for b in range(2):
            d_means[b, a] = means[a+1, b] - means[0, b]
            d_std[b, a] = stand[a+1, b] - stand[0, b]
    d_means = np.reshape(d_means, (2, 9, 34))
    d_std = np.reshape(d_std, (2, 9, 34))
    # change sign of negative latitudes
    d_means[0, :] = d_means[0, :] * -1

    return daily, means, d_means, stand, d_std


def main2KHadstrength():
    def getHadST(strfcn, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            hadstr = np.min(np.mean(strfcn[22:31, 33:], axis=0))
        else:
            hadstr = np.max(np.mean(strfcn[22:31, :], axis=0))
        return hadstr

    daily = np.zeros((307, 72, 2))
    # run, day, winter/summer
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        hadext = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadST(strfcn[t, :], hem='NH')
            hadext[t, 1] = getHadST(strfcn[t, :], hem='SH')

        daily[i, :, 1], daily[i, :, 0] = hadext[:, 0], hadext[:, 1]
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    means = I2Kfull(means)
    d_means = np.zeros((2, 306))
    for a in range(306):
        for b in range(2):
            d_means[b, a] = means[a+1, b] - means[0, b]

    d_means = np.reshape(d_means, (2, 9, 34))
    d_means[1, :] = d_means[1, :] * -1

    return daily, means, d_means


def main2Kitczwidth():
    def getHadEx(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(50.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
            lat1 = lat[idx_min:idx_max+1]

        a = np.argmin(strfcn1)           # Find the lat index of global minimum

        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))
    def getHadST(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            hadstr = lat[np.argmin(np.mean(strfcn[22:31, 33:], axis=0)) + 33]
        else:
            hadstr = lat[np.argmax(np.mean(strfcn[22:31, :], axis=0))]
        return hadstr

    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')

    daily = np.zeros((307, 72))
    # run, day, winter/summer
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        hadext = np.zeros((72, 2))
        hadext1 = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadST(strfcn[t, :], lat, hem='NH')
            hadext[t, 1] = getHadST(strfcn[t, :], lat, hem='SH')
            hadext1[t, 0] = getHadEx(strfcn[t, :], lat, hem='NH')
            hadext1[t, 1] = getHadEx(strfcn[t, :], lat, hem='SH')

        daily[i, :] = (hadext[:, 0] - hadext[:, 1]) / (hadext1[:, 0] - hadext1[:, 1])
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    means = I2Kfull(means)
    d_means = np.zeros((306))
    for a in range(306):
        d_means[a] = means[a+1] - means[0]

    d_means = np.reshape(d_means, (9, 34))
    return daily, means, d_means


def main2KITCZ_vertav():
    def getHadEx(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(50.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            lat1 = lat[idx_min:idx_max+1]
        elif hem == 'ITCZ':
            idx_min = getIdx(0, lat)
            idx_max = getIdx(25.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
            lat1 = lat[idx_min:idx_max+1]

        a = np.argmin(strfcn1)           # Find the lat index of global minimum

        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'ITCZ':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))

    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')

    daily = np.zeros((307, 72, 3))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        hadext = np.zeros((72, 3))
        for t in range(72):
            hadext[t, 0] = getHadEx(strfcn[t, :], lat, hem='NH')
            hadext[t, 1] = getHadEx(strfcn[t, :], lat, hem='ITCZ')
            hadext[t, 2] = getHadEx(strfcn[t, :], lat, hem='SH')

        daily[i, :, 2], daily[i, :, 1], daily[i, :, 0] = (hadext[:, 0]-hadext[:, 1]), hadext[:, 1], (hadext[:, 1]-hadext[:, 2])
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    stand = np.std(daily, axis=1)
    means = I2Kfull(means)
    stand = I2Kfull(stand)
    d_means = np.zeros((3, 306))
    d_std = np.zeros((3, 306))
    for a in range(306):
        for b in range(3):
            d_means[b, a] = means[a+1, b] - means[0, b]
            d_std[b, a] = stand[a+1, b] - stand[0, b]
    d_means = np.reshape(d_means, (3, 9, 34))
    d_std = np.reshape(d_std, (3, 9, 34))
    # change sign of negative latitudes
    # d_means[0, :] = d_means[0, :] * -1

    return daily, means, d_means, stand, d_std


def hadextentpredict():

    # Try it using a constant diffusivity.
    # If it looks a total mess, then diagnose the diffusivity from the GCM
    # using D = [v'T']/ [grad(T)] where [] is a near-surface
    # and midlatitude (35N-65N) average.

    def ebm_strklat2(G, delta, hadEx, D):
        # Temperature gradient at the Hadley Terminus: G
        # The pole to equator thermal contrast in radiative equilibrium: delta
        # Hadley Terminus in degrees: hadEx
        # Hugh - think delta=95K in my experiments
        # D = 2.1e6  # Diffusivity
        R = 6.365e6  # Radius of Earth
        A = 1/(50*86400)  # Inverse radiative timescale
        yh = (abs(hadEx) * np.pi/180 * R)  # Hadley Terminus in meters north
        alpha = 1/(R*np.pi/2 - yh)
        zeta = np.pi**2 * alpha**2 * D/A
        rad = (1/8)*delta*R*(np.pi*R - 2*yh)*(np.sin(2*yh/R)/(np.pi*R*yh - yh**2) - 32*np.sin(np.pi*10/180)*np.cos(yh/R)/((3*np.pi*R - 2*yh)*(np.pi*R + 2*yh)))
        psi = 3/4 * (1 + zeta)/(4 + zeta)/((3 * np.pi**2)/(alpha) * (rad)/(4+zeta) * 1/(abs(G/alpha**2)) - 1)
        strklat = (yh + 0.5*(np.sqrt(psi**2 + 8/(np.pi**2)) - psi)/alpha)/R*180/np.pi
        return strklat

    def getHadEx(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(50.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
            lat1 = lat[idx_min:idx_max+1]

        a = np.argmin(strfcn1)           # Find the lat index of global minimum

        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))

    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')
    R = 6.365e6
    daily = np.zeros((307, 72, 2))
    dailyG = np.zeros((307, 72, 2))
    dailyD = np.zeros((307, 72, 2))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        temp = indiv('2Krun', str(i), 'temp')[12:, :, :]
        vT = indiv('2Krun', str(i), 'temp_eddy_mrdnl_flux')[12:, :, :]
        hadext = np.zeros((72, 2))
        t_grad = np.zeros((72, 2))
        diff = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadEx(strfcn[t, :], lat, hem='SH')
            hadext[t, 1] = getHadEx(strfcn[t, :], lat, hem='NH')
            # code to get T gradient
            a0 = getIdx(hadext[t, 0], lat)
            a1 = getIdx(hadext[t, 1], lat)
            t_mean = np.mean(temp[t, 22:31, :], axis=0)
            t_mean1 = np.mean(temp[t, 22:31, :], axis=0)
            t_grad[t, 0] = (t_mean[a0-2] - t_mean[a0+3])/((lat[a0-2]-lat[a0+3])*R*np.pi/180)
            t_grad[t, 1] = (t_mean[a1-2] - t_mean[a1+3])/((lat[a1-2]-lat[a1+3])*R*np.pi/180)
            diff[t, 0] = np.mean(vT[t, 28:31, 8:a0])/((t_mean1[8] - t_mean1[a0])/((lat[8]-lat[a0])*R*np.pi/180))
            diff[t, 1] = np.mean(vT[t, 28:31, a1:56])/((t_mean1[a1] - t_mean1[56])/((lat[a1]-lat[56])*R*np.pi/180))
        dailyG[i, :, 0], dailyG[i, :, 1] = t_grad[:, 0], t_grad[:, 1]
        dailyD[i, :, 0], dailyD[i, :, 1] = abs(diff[:, 0]), abs(diff[:, 1])
        daily[i, :, 0], daily[i, :, 1] = hadext[:, 0], hadext[:, 1]
        print('Completed run' + str(i))

    means = np.mean(daily, axis=1)
    means = I2Kfull(means)
    G = np.mean(dailyG, axis=1)
    G = I2Kfull(G)
    D = np.mean(dailyD, axis=1)
    D = I2Kfull(D)

    strklats = np.zeros((2, 307))
    for i in range(307):
        for a in range(2):
                strklats[a, i] = ebm_strklat2(G[i, a], 95, means[i, a], D[i, a])
    strklats[0, :] = strklats[0, :]*-1
    return strklats, means, G


def main2KHad_vertav_jetdiff():
    def getHadEx(strfcn, lat, hem='NH'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(50.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
            lat1 = lat[idx_min:idx_max+1]

        a = np.argmin(strfcn1)           # Find the lat index of global minimum

        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))

    sigma = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                   '2Krun/run0/history/day3630h00.nc', 'sigma')
    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')
    gridlat = np.concatenate((lat[0:32:2], lat[31:33], lat[33::2]))
    gridsigma = np.zeros((9))
    gridsigma[0], gridsigma[2], gridsigma[4] = sigma[36], sigma[30], sigma[26]
    gridsigma[6], gridsigma[8] = sigma[22], sigma[16]
    gridsigma[1], gridsigma[3], gridsigma[5] = sigma[33], sigma[28], sigma[24]
    gridsigma[7] = sigma[19]

    daily = np.zeros((307, 72, 2))
    daily_normj = np.zeros((307, 72, 2))
    daily_normh = np.zeros((307, 72, 2))
    # u850 = np.zeros((307, 2160, 64))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        hadext = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadEx(strfcn[t, :], lat, hem='NH')
            hadext[t, 1] = getHadEx(strfcn[t, :], lat, hem='SH')
        udaily1 = indiv('2Krun', str(i), 'u')[12:, 31, :]
        js, jl = jetindicesdaily(udaily1, lat)
        daily[i, :, 1], daily[i, :, 0] = jl[:, 1] - hadext[:, 0], jl[:, 0] - hadext[:, 1]
        daily_normj[i, :, 1], daily_normj[i, :, 0] = jl[:, 1], jl[:, 0]
        daily_normh[i, :, 1], daily_normh[i, :, 0] = hadext[:, 0], hadext[:, 1]
        print('Completed polyfit for run' + str(i))

    means = np.mean(daily, axis=1)
    meansj= np.mean(daily_normj, axis=1)
    meansh = np.mean(daily_normh, axis=1)
    sd = np.std(daily, axis=1)
    means = I2Kfull(means)
    meansj = I2Kfull(meansj)
    meansh = I2Kfull(meansh)
    sd = I2Kfull(sd)
    d_means = np.zeros((2, 306))
    d_meansj = np.zeros((2, 306))
    d_meansh = np.zeros((2, 306))
    d_sd = np.zeros((2, 306))
    for a in range(306):
        for b in range(2):
                d_means[b, a] = means[a+1, b] - means[0, b]
                d_meansj[b, a] = meansj[a+1, b] - meansh[0, b]
                d_meansh[b, a] = meansh[a+1, b] - meansj[0, b]
                d_sd[b, a] = sd[a+1, b] - sd[0, b]

    d_means = np.reshape(d_means, (2, 9, 34))
    d_meansj = np.reshape(d_meansj, (2, 9, 34))
    d_meansh = np.reshape(d_meansh, (2, 9, 34))
    d_sd = np.reshape(d_sd, (2, 9, 34))

    # change sign of negative latitudes
    d_means[0, :, :] = d_means[0, :, :] * -1
    d_meansj[0, :, :] = d_meansj[0, :, :] * -1
    d_meansh[0, :, :] = d_meansh[0, :, :] * -1

    d_meansj = d_means / d_meansj
    d_meansh = d_means / d_meansh
    return d_means, d_sd, d_meansj, d_meansh


def getFerEx(strfcn, lat, hem='NH'):
    # Isolate the northern hemisphere
    if hem == 'NH':
        idx_min = getIdx(20.0, lat)
        idx_max = getIdx(90.0, lat)
        strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
        lat1 = lat[idx_min:idx_max+1]
    else:
        idx_max = getIdx(-0.5, lat)
        idx_min = getIdx(-90.0, lat)
        strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
        lat1 = lat[idx_min:idx_max+1]

    # a = np.argmin(strfcn1)           # Find the lat index of global minimum

    # Change of sign of the streamfunction
    if hem == 'NH':
        a = -1
        aa = a - 1
        while np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
            aa = aa - 1
        hadEx = np.interp(0.0, np.asarray([strfcn1[aa],
                                           strfcn1[aa+1]]),
                          np.asarray([lat1[aa], lat1[aa+1]]))
    elif hem == 'SH':
        a = 0
        aa = a + 1
        while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
            aa = aa + 1
        hadEx = np.interp(0.0, np.asarray([strfcn1[aa],
                                           strfcn1[aa-1]]),
                          np.asarray([lat1[aa], lat1[aa-1]]))
    return hadEx


def getHadEx(strfcn, lat, hem='NH'):
    # Isolate the northern hemisphere
    if hem == 'NH':
        idx_min = getIdx(20.0, lat)
        idx_max = getIdx(50.0, lat)
        strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
        lat1 = lat[idx_min:idx_max+1]
    else:
        idx_max = getIdx(-0.5, lat)
        idx_min = getIdx(-90.0, lat)
        strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
        lat1 = lat[idx_min:idx_max+1]

    a = np.argmin(strfcn1)           # Find the lat index of global minimum

    # Change of sign of the streamfunction
    if hem == 'NH':
        aa = a - 1
        while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
            aa = aa - 1
        hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                           strfcn1[aa]]),
                          np.asarray([lat1[aa+1], lat1[aa]]))
    elif hem == 'SH':
        aa = a + 1
        while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
            aa = aa + 1
        hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                           strfcn1[aa]]),
                          np.asarray([lat1[aa-1], lat1[aa]]))
    return hadEx


def getIdx(x0, x):
    y = x.tolist()
    return y.index(min(y, key=lambda y: abs(y-x0)))


def getTvT():
    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')
    dailyT = np.zeros((307, 72, 2))
    dailyvT = np.zeros((307, 72, 2))
    dailyT_f = np.zeros((307, 72, 2))
    dailyvT_f = np.zeros((307, 72, 2))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        temp = indiv('2Krun', str(i), 'temp')[12:, :, :]
        vTs = indiv('2Krun', str(i), 'temp_eddy_mrdnl_flux')[12:, :, :]
        hadext = np.zeros((72, 2))
        ferext = np.zeros((72, 2))
        t_grad = np.zeros((72, 2))
        vT = np.zeros((72, 2))
        t_grad_f = np.zeros((72, 2))
        vT_f = np.zeros((72, 2))
        for t in range(72):
            hadext[t, 0] = getHadEx(strfcn[t, :], lat, hem='SH')
            hadext[t, 1] = getHadEx(strfcn[t, :], lat, hem='NH')
            ferext[t, 0] = getFerEx(strfcn[t, :], lat, hem='SH')
            ferext[t, 1] = getFerEx(strfcn[t, :], lat, hem='NH')
            # code to get T gradient
            a0 = getIdx(hadext[t, 0], lat)
            a1 = getIdx(hadext[t, 1], lat)
            a2 = getIdx(ferext[t, 0], lat)
            a3 = getIdx(ferext[t, 1], lat)
            t_mean = np.mean(temp[t, 22:31, :], axis=0)
            t_grad[t, 0] = (t_mean[a0-2] - t_mean[a0+3])/((lat[a0-2]-lat[a0+3]))
            t_grad[t, 1] = (t_mean[a1-2] - t_mean[a1+3])/((lat[a1-2]-lat[a1+3]))
            vT[t, 0] = np.mean(vTs[t, 22:31, a0], axis=0)
            vT[t, 1] = np.mean(vTs[t, 22:31, a1], axis=0)
            t_grad_f[t, 0] = (t_mean[a2-1] - t_mean[a2+2])/((lat[a2-1]-lat[a2+2]))
            if a3 < 61:
                t_grad_f[t, 1] = (t_mean[a3-1] - t_mean[a3+2])/((lat[a3-1]-lat[a3+2]))
            else:
                t_grad_f[t, 1] = (t_mean[a3-2] - t_mean[a3])/((lat[a3-2]-lat[a3]))
            vT_f[t, 0] = np.mean(vTs[t, 22:31, a2], axis=0)
            vT_f[t, 1] = np.mean(vTs[t, 22:31, a3], axis=0)
        dailyT[i, :, 0], dailyT[i, :, 1] = t_grad[:, 0], t_grad[:, 1]
        dailyvT[i, :, 0], dailyvT[i, :, 1] = abs(vT[:, 0]), abs(vT[:, 1])
        dailyT_f[i, :, 0], dailyT_f[i, :, 1] = t_grad_f[:, 0], t_grad_f[:, 1]
        dailyvT_f[i, :, 0], dailyvT_f[i, :, 1] = abs(vT_f[:, 0]), abs(vT_f[:, 1])
        print('Completed run' + str(i))

    T = np.mean(dailyT, axis=1)
    T = I2Kfull(T)
    H = np.mean(dailyvT, axis=1)
    H = I2Kfull(H)
    Tf = np.mean(dailyT_f, axis=1)
    Tf = I2Kfull(Tf)
    Hf = np.mean(dailyvT_f, axis=1)
    Hf = I2Kfull(Hf)
    return T, Tf, H, Hf


def main2Ksub_jet():
    def getHadEx(strfcn, lat, hem='NH', method='strm'):
        # Isolate the northern hemisphere
        if hem == 'NH':
            idx_min = getIdx(20.0, lat)
            idx_max = getIdx(60.0, lat)
            if method == 'strm':
                strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0) * -1.0
            else:
                strfcn1 = strfcn[16, idx_min:idx_max+1]
            lat1 = lat[idx_min:idx_max+1]
        else:
            idx_max = getIdx(-0.5, lat)
            idx_min = getIdx(-90.0, lat)
            if method == 'strm':
                strfcn1 = np.mean(strfcn[22:31, idx_min:idx_max+1], axis=0)
            else:
                strfcn1 = strfcn[16, idx_min:idx_max+1]
            lat1 = lat[idx_min:idx_max+1]

        a = np.argmin(strfcn1)             # Find the lat index of global minimum

        # Change of sign of the streamfunction
        if hem == 'NH':
            aa = a - 1
            while aa > 0 and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa - 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa+1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa+1], lat1[aa]]))
        elif hem == 'SH':
            aa = a + 1
            while aa < strfcn1.size and np.sign(strfcn1[a]) == np.sign(strfcn1[aa]):
                aa = aa + 1
            hadEx = np.interp(0.0, np.asarray([strfcn1[aa-1],
                                               strfcn1[aa]]),
                              np.asarray([lat1[aa-1], lat1[aa]]))
        return hadEx

    def getIdx(x0, x):
        y = x.tolist()
        return y.index(min(y, key=lambda y: abs(y-x0)))

    lat = ncread('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                 '2Krun/run0/history/day3630h00.nc', 'lat')

    def inter(u, lat):
        x0 = int(np.floor(lat))
        x1 = int(np.ceil(lat))
        y0 = u[:, x0]
        y1 = u[:, x1]
        y = y0 + (lat-x0)*(y1-y0)/(x1-x0)
        return y

    daily_shear = np.zeros((307, 72, 2))
    daily_umax = np.zeros((307, 72, 2))
    daily_duvmax = np.zeros((307, 72, 2))
    # run, day, winter/summer, speed/lat
    for i in range(307):
        strfcn = indiv('2Krun', str(i), 'streamfctn')[12:, :, :]
        u = indiv('2Krun', str(i), 'u')[12:, :, :]
        uv = indiv('2Krun', str(i), 'u_eddy_mrdnl_flux')[12:, :, :]
        duv = np.gradient(uv, axis=2)
        for t in range(72):
            daily_umax[i, t, 1] = np.max(inter(u[t, :25], getHadEx(strfcn[t, :], lat, hem='NH')))
            daily_umax[i, t, 0] = np.max(inter(u[t, :25], getHadEx(strfcn[t, :], lat, hem='SH')))
            daily_shear[i, t, 1] = np.max(u[t, :25, 32:]-u[t, 31, 32:])
            daily_shear[i, t, 0] = np.max(u[t, :25, :32]-u[t, 31, :32])
            daily_duvmax[i, t, 1] = np.max(inter(u[t, :25], getHadEx(duv[t, :], lat, hem='NH', method='duv')))
            daily_duvmax[i, t, 0] = np.max(inter(u[t, :25], getHadEx(duv[t, :], lat, hem='SH', method='duv')))
            '''
            daily_umax[i, t, 1] = np.max(inter(u[t], getHadEx(strfcn[t, :], lat, hem='NH')))
            daily_umax[i, t, 0] = np.max(inter(u[t], getHadEx(strfcn[t, :], lat, hem='SH')))
            daily_shear[i, t, 1] = lat[np.remainder(np.argmax(u[t, :25, 32:]-u[t, 31, 32:]), 32)+32]
            daily_shear[i, t, 0] = lat[np.remainder(np.argmax(u[t, :25, :32]-u[t, 31, :32]), 32)]
            daily_duvmax[i, t, 1] = getHadEx(duv[t, :], lat, hem='NH', method='duv')
            daily_duvmax[i, t, 0] = getHadEx(duv[t, :], lat, hem='SH', method='duv')
            '''
        print('Completed run' + str(i))

    meansu = np.mean(daily_shear, axis=1)
    meansu = I2Kfull(meansu)
    d_meansu = np.zeros((2, 306))
    meansum = np.mean(daily_umax, axis=1)
    meansum = I2Kfull(meansum)
    d_meansum = np.zeros((2, 306))
    meansduv = np.mean(daily_duvmax, axis=1)
    meansduv = I2Kfull(meansduv)
    d_meansduv = np.zeros((2, 306))
    #  d_sd = np.zeros((2, 306))
    for a in range(306):
        for b in range(2):
            d_meansu[b, a] = meansu[a+1, b] - meansu[0, b]
            d_meansum[b, a] = meansum[a+1, b] - meansum[0, b]
            d_meansduv[b, a] = meansduv[a+1, b] - meansduv[0, b]

    d_meansu = np.reshape(d_meansu, (2, 9, 34))
    d_meansum = np.reshape(d_meansum, (2, 9, 34))
    d_meansduv = np.reshape(d_meansduv, (2, 9, 34))
    # change sign of negative latitudes
    d_meansu[0, :] = d_meansu[0, :] * -1
    d_meansum[0, :] = d_meansum[0, :] * -1
    d_meansduv[0, :] = d_meansduv[0, :] * -1

    return daily_shear, d_meansu, daily_duvmax, d_meansduv, daily_umax, d_meansum
