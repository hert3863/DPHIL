#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 16:00:29 2018

@author: bakerh
"""


import numpy as np
import scipy as sp
import scipy.stats as stats
import matplotlib.pyplot as plt


def plot_ci_manual(t, s_err, n, x, x2, y2, ax=None):
    """Return an axes of confidence bands using a simple approach.

    Notes
    -----
    .. math:: \left| \: \hat{\mu}_{y|x0} - \mu_{y|x0} \: \right| \; \leq \; T_{n-2}^{.975} \; \hat{\sigma} \; \sqrt{\frac{1}{n}+\frac{(x_0-\bar{x})^2}{\sum_{i=1}^n{(x_i-\bar{x})^2}}}
    .. math:: \hat{\sigma} = \sqrt{\sum_{i=1}^n{\frac{(y_i-\hat{y})^2}{n-2}}}

    References
    ----------
    .. [1]: M. Duarte.  "Curve fitting," JUpyter Notebook.
       http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/CurveFitting.ipynb

    """
    if ax is None:
        ax = plt.gca()

    ci = t*s_err*np.sqrt(1/n + (x2-np.mean(x))**2/np.sum((x-np.mean(x))**2))
    ax.fill_between(x2, y2+ci, y2-ci, color="#b9cfe7", edgecolor="")

    return ax


def regplotter(x, y):
    # Modeling with Numpy
    p, cov = np.polyfit(x, y, 1, cov=True)        # parameters and covariance from of the fit
    y_model = np.polyval(p, x)                    # model using the fit parameters; NOTE: parameters here are coefficients

    # Statistics
    n = y.size                              # number of observations
    m = p.size                                    # number of parameters
    DF = n - m                                    # degrees of freedom
    t = stats.t.ppf(5/6, n - m)                  # used for CI and PI bands

    # Estimates of Error in Data/Model
    resid = y - y_model
    chi2 = np.sum((resid/y_model)**2)             # chi-squared; estimates error in data
    chi2_red = chi2/(DF)                          # reduced chi-squared; measures goodness of fit
    s_err = np.sqrt(np.sum(resid**2)/(DF))        # standard deviation of the error

    # Plotting ---------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 6))

    # Data
    ax.plot(x, y, "o", color="#b9cfe7", markersize=5,
            markeredgewidth=1, markeredgecolor="b", markerfacecolor="b")

    # Fit
    ax.plot(x, y_model, "-", color="0.1", linewidth=1.5, alpha=0.5,
            label="Fit")

    x2 = np.linspace(np.min(x), np.max(x), 100)
    y2 = np.linspace(np.min(y_model), np.max(y_model), 100)

    # Confidence Interval (select one)
    plot_ci_manual(t, s_err, n, x, x2, y2, ax=ax)
    # plot_ci_bootstrap(n, x, y, resid, ax=ax)

    # Prediction Interval
    pi = t*s_err*np.sqrt(1+1/n+(x2-np.mean(x))**2/np.sum((x-np.mean(x))**2))
    ax.fill_between(x2, y2+pi, y2-pi, color="None", linestyle="--")
    ax.plot(x2, y2-pi, "--", color="0.5", label="66% Prediction Limits")
    ax.plot(x2, y2+pi, "--", color="0.5")

    # Figure Modifications ---------------------------------------------------
    # Borders
    ax.spines["top"].set_color("0.5")
    ax.spines["bottom"].set_color("0.5")
    ax.spines["left"].set_color("0.5")
    ax.spines["right"].set_color("0.5")
    ax.get_xaxis().set_tick_params(direction="out")
    ax.get_yaxis().set_tick_params(direction="out")
    ax.xaxis.tick_bottom()
    ax.yaxis.tick_left()

    # Labels
    #plt.title("Fit Plot for Weight", fontsize="14", fontweight="bold")
    plt.xlabel("Land-sea temperature contrast change (K)")
    plt.ylabel("ECS (K)")
    plt.xlim(np.min(x)-.1, np.max(x)+.1)

    # Custom legend
    handles, labels = ax.get_legend_handles_labels()
    display = (0, 1)
    anyArtist = plt.Line2D((0, 1), (0, 0), color="#b9cfe7")  # Create custom artists
    legend = plt.legend(
              [handle for i, handle in enumerate(handles) if i in display]+[anyArtist],
              [label for i, label in enumerate(labels) if i in display]+["66% Confidence Limits"],
              loc=9, bbox_to_anchor=(0, -0.21, 1., .102), ncol=3, mode="expand")  
    frame = legend.get_frame().set_edgecolor("0.5")

    # Save Figure
    #plt.tight_layout()
    # plt.savefig("filename.png", bbox_extra_artists=(legend,), bbox_inches="tight")

    plt.show()


def cmiptas():
    from scipy import interpolate
    from netcdfread import ncread
    import glob
    import numpy as np

    def lscalc(t_w, lsm1, meshlatweight, r1=56, r2=89):
        lsm1 = lsm1[r1:r2, :]
        t_w = t_w[r1:r2, :]
        meshlatweight = meshlatweight[r1:r2, :]
        t_l = np.mean(t_w[lsm1 > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
        t_o = np.mean(t_w[lsm1 < 0.5]) / np.mean(meshlatweight[lsm1 < 0.5])
        return t_l, t_o

    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    lat_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'latitude')
    lon_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'longitude')
    meshlat = np.zeros([np.ma.size(lon_n96), np.ma.size(lat_n96)])
    meshlat[:, :] = lat_n96
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    models = ['ACCESS1-0_', 'ACCESS1-3_', 'bcc-csm1-1_', 'bcc-csm1-1-m_',
              'BNU-ESM_', 'CanESM2_', 'CCSM4_', 'CNRM-CM5_', 'CSIRO-Mk3-6-0_',
              'EC-EARTH_', 'FGOALS-g2_', 'FGOALS-s2_', 'GFDL-CM3_',
              'GFDL-ESM2G_', 'GFDL-ESM2M_', 'GISS-E2-H_', 'GISS-E2-R_',
              'HadGEM2-ES', 'inmcm4_', 'IPSL-CM5A-LR_', 'IPSL-CM5A-MR_',
              'IPSL-CM5B-LR_', 'MPI-ESM-LR_', 'MPI-ESM-MR_', 'MPI-ESM-P_',
              'NorESM1-M_']
    dT = np.zeros(26)
    aT = np.zeros(26)
    aTt = np.zeros(26)
    for j, model in enumerate(models):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/abrupt4xco2/tas_Amon_' + model + '*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/piControl/tas_Amon_' + model + '*')
        data4 = ncread(a[0], 'tas')
        datac = ncread(b[0], 'tas')
        if len(a) > 1:
            for i in range(len(a)-1):
                data4 = np.concatenate((data4, ncread(a[i+1], 'tas')), axis=0)
        if len(b) > 1:
            for i in range(len(b)-1):
                datac = np.concatenate((datac, ncread(b[i+1], 'tas')), axis=0)
        data4 = np.squeeze(np.mean(data4[1200:1800], axis=0))
        datac = np.squeeze(np.mean(datac[1200:], axis=0))
        lat = ncread(a[0], 'lat')
        lon = ncread(a[0], 'lon')
        f = interpolate.interp2d(lon, lat[::-1], data4)
        g = interpolate.interp2d(lon, lat[::-1], datac)
        data4_n96_w = f(lon_n96, lat_n96)*meshlatweight
        datac_n96_w = g(lon_n96, lat_n96)*meshlatweight
        t4 = np.mean(data4_n96_w)/np.mean(meshlatweight)
        tc = np.mean(datac_n96_w)/np.mean(meshlatweight)
        dt = t4 - tc
        t4_l, t4_o = lscalc(data4_n96_w, lsm1, meshlatweight, 0, 145)
        tc_l, tc_o = lscalc(datac_n96_w, lsm1, meshlatweight, 0, 145)
        t4t_l, t4t_o = lscalc(data4_n96_w, lsm1, meshlatweight, 56, 89)
        tct_l, tct_o = lscalc(datac_n96_w, lsm1, meshlatweight, 56, 89)
        ta = (t4_l-tc_l) - (t4_o-tc_o)
        tat = (t4t_l-tct_l) - (t4t_o-tct_o)
        dT[j] = dt
        aT[j] = ta
        aTt[j] = tat
        print(model)
    return dT, aT, aTt


def mean_4xco2(var='va'):
    from scipy import interpolate
    from netcdfread import ncread
    import glob
    import numpy as np
    from btmodel import maplot_ll
    lat_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'latitude')
    lon_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'longitude')
    a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/AMIP/' + var + '/*'))
    b = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/AMIP4xCO2/' + var + '/*'))
    sindices = np.zeros((30*3))
    for i in range(30):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    sindices = sindices.astype(int)
    v_means = np.zeros((9, 145, 192))
    for i in range(9):
        st = 0
        if i == 3:
            st = 4
        if i == 1:
            st = 348
        if var == 'pr':
            v_c = ncread(a[i], var)[st+sindices]
            v_f = ncread(b[i], var)[st+sindices]
        else:
            v_c = ncread(a[i], var)[st+sindices, 9]
            v_f = ncread(b[i], var)[st+sindices, 9]
        dv = np.mean(v_f-v_c, axis=0)
        lat = ncread(a[i], 'lat')
        lon = ncread(b[i], 'lon')
        maplot_ll(dv, lat, lon, 3)
        f = interpolate.interp2d(lon, lat[::-1], dv)
        v_res = f(lon_n96, lat_n96)
        v_means[i] = v_res
    return v_means


def mean_4K(var='va'):
    from scipy import interpolate
    from netcdfread import ncread
    import glob
    import numpy as np
    lat_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'latitude')
    lon_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'longitude')
    a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/AMIP/' + var + '/*'))
    b = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/AMIP4K/' + var + '/*'))
    sindices = np.zeros((30*3))
    for i in range(30):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    sindices = sindices.astype(int)
    v_means = np.zeros((9, 145, 192))
    for i in range(9):
        st = 0
        if i == 3:
            st = 4
        if i == 1:
            st = 348
        if var == 'pr':
            v_c = ncread(a[i], var)[st+sindices]
            v_f = ncread(b[i], var)[st+sindices]
        else:
            v_c = ncread(a[i], var)[st+sindices, 9]
            v_f = ncread(b[i], var)[st+sindices, 9]
        dv = np.mean(v_f-v_c, axis=0)
        lat = ncread(a[i], 'lat')
        lon = ncread(b[i], 'lon')

        f = interpolate.interp2d(lon, lat[::-1], dv)
        v_res = f(lon_n96, lat_n96)
        v_means[i] = v_res
    return v_means


def lsc(tbar, r1, r2):
    from netcdfread import ncread
    import numpy as np
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    t_nat = np.mean(tbar['All-Nat'], axis=0)
    t_sst = np.mean(tbar['GHG-Nat'], axis=0)
    t_ghg = np.mean(tbar['SST-Nat'], axis=0)
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))

    t_nat_w = t_nat * meshlatweight
    t_sst_w = t_sst * meshlatweight
    t_ghg_w = t_ghg * meshlatweight

    def lscalc(t_w, lsm1, meshlatweight, r1=56, r2=89):
        lsm1 = lsm1[r1:r2, :]
        t_w = t_w[r1:r2, :]
        meshlatweight = meshlatweight[r1:r2, :]
        t_l = np.mean(t_w[lsm1 > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
        t_o = np.mean(t_w[lsm1 < 0.5]) / np.mean(meshlatweight[lsm1 < 0.5])
        return t_l, t_o

    t_nat_l, t_nat_o = lscalc(t_nat_w, lsm1, meshlatweight, r1, r2)
    t_sst_l, t_sst_o = lscalc(t_sst_w, lsm1, meshlatweight, r1, r2)
    t_ghg_l, t_ghg_o = lscalc(t_ghg_w, lsm1, meshlatweight, r1, r2)
    return t_nat_l, t_nat_o, t_sst_l, t_sst_o, t_ghg_l, t_ghg_o


def cmipw5(v_had):
    import glob
    from netcdfread import ncread
    from scipy import interpolate
    lat_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'latitude')
    lon_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'longitude')

    models = ['ACCESS1-0_', 'ACCESS1-3_', 'bcc-csm1-1_', 'bcc-csm1-1-m_',
              'BNU-ESM_', 'CanESM2_', 'CCSM4_', 'CESM1-BGC_', 'CESM1-CAM5_',
              'CMCC-CESM_', 'CMCC-CM_', 'CMCC-CMS_', 'CNRM-CM5_',
              'CSIRO-Mk3-6-0_',
              'EC-EARTH_', 'FGOALS-g2_', 'FIO-ESM_', 'GFDL-CM3_',
              'GFDL-ESM2G_', 'GFDL-ESM2M_', 'GISS-E2-H-CC_', 'GISS-E2-H_',
              'GISS-E2-R-CC_', 'GISS-E2-R_', 'HadGEM2-AO_', 'HadGEM2-CC_',
              'HadGEM2-ES_', 'inmcm4_', 'IPSL-CM5A-LR_', 'IPSL-CM5A-MR_',
              'IPSL-CM5B-LR_', 'MPI-ESM-LR_', 'MPI-ESM-MR_',
              'NorESM1-ME_', 'NorESM1-M_']
    indices = {}
    for k, model in enumerate(models):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/historical/va/va_Amon_' + model + '*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/rcp85/va/va_Amon_' + model + '*')
        m = 6
        e = 156
        time = np.arange(1850, 2101)
        if model == 'bcc-csm1-1_':
            time = np.arange(1850, 2301)
        elif model == 'IPSL-CM5A-LR_':
            time = np.arange(1850, 2301)
        elif model == 'CCSM4_':
            time = np.arange(1850, 2301)
        elif model == 'CNRM-CM5_':
            time = np.arange(1850, 2301)
        elif model == 'CSIRO-Mk3-6-0_':
            time = np.arange(1850, 2301)
        elif model == 'GISS-E2-R_':
            time = np.arange(1850, 2301)
        elif model == 'GISS-E2-H_':
            time = np.arange(1850, 2301)
        elif model == 'MPI-ESM-LR_':
            time = np.arange(1850, 2301)
        elif model == 'GFDL-CM3_':
            e = 146
            time = np.arange(1860, 2101)
        elif model == 'GFDL-ESM2G_':
            e = 145
            time = np.arange(1861, 2101)
        elif model == 'GFDL-ESM2M_':
            e = 145
            time = np.arange(1861, 2101)
        elif model == 'HadGEM2-CC_':
            m = 5
            e = 146
            time = np.arange(1860, 2100)
        elif model == 'HadGEM2-AO_':
            e = 146
            time = np.arange(1860, 2100)
        elif model == 'HadGEM2-ES_':
            m = 5
            e = 146
            time = np.arange(1860, 2300)
        elif model == 'FGOALS-g2_':
            m = 6
            e = 106
            time = np.arange(1900, 2102)

        va = ncread(a[0], 'va')[m::12, 9]
        va1 = ncread(b[0], 'va')[m::12, 9]
        if len(a) > 1:
            for i in range(len(a)-1):
                va = np.concatenate((va, ncread(a[i+1], 'va')[m::12, 9]), axis=0)
        if len(b) > 1:
            for i in range(len(b)-1):
                va1 = np.concatenate((va1, ncread(b[i+1], 'va')[m::12, 9]), axis=0)
        lat = ncread(a[0], 'lat')
        lon = ncread(a[0], 'lon')
        va = va[:e]
        va_j2 = np.concatenate((va, va1), axis=0)
        va_j2_96 = np.zeros((len(va_j2), 145, 192))
        for j in range(len(va_j2)):
            f = interpolate.interp2d(lon, lat[::-1], va_j2[j])
            va_j2_96[j] = f(lon_n96, lat_n96)
        meshlat = np.zeros([np.ma.size(lon_n96), np.ma.size(lat_n96)])
        meshlat[:, :] = lat_n96
        meshlatweight = np.cos(meshlat.transpose() * np.pi/180)
        index_j2 = np.zeros(len(va_j2))
        for i in range(len(va_j2)):
            index_j2[i] = np.sum(v_had[12:49]*va_j2_96[i, 12:49]*meshlatweight[12:49])/(np.sqrt(np.sum(v_had[12:49]*v_had[12:49]*meshlatweight[12:49]))*np.sqrt(np.sum(va_j2_96[i, 12:49]*va_j2_96[i, 12:49]*meshlatweight[12:49])))
        indices[model] = time, index_j2
        print(model)
    return indices


def plotcmip5(indices):
    import matplotlib.pyplot as plt
    models = ['ACCESS1-0_', 'ACCESS1-3_', 'bcc-csm1-1_', 'bcc-csm1-1-m_',
          'BNU-ESM_', 'CanESM2_', 'CCSM4_', 'CESM1-BGC_', 'CESM1-CAM5_',
          'CMCC-CESM_', 'CMCC-CM_', 'CMCC-CMS_', 'CNRM-CM5_',
          'CSIRO-Mk3-6-0_',
          'EC-EARTH_', 'FGOALS-g2_', 'FIO-ESM_', 'GFDL-CM3_',
          'GFDL-ESM2G_', 'GFDL-ESM2M_', 'GISS-E2-H-CC_', 'GISS-E2-H_',
          'GISS-E2-R-CC_', 'GISS-E2-R_', 'HadGEM2-AO_', 'HadGEM2-CC_',
          'HadGEM2-ES_', 'inmcm4_', 'IPSL-CM5A-LR_', 'IPSL-CM5A-MR_',
          'IPSL-CM5B-LR_', 'MPI-ESM-LR_', 'MPI-ESM-MR_',
          'NorESM1-ME_', 'NorESM1-M_']

    for model in models:
        if model != 'EC-EARTH_':
            plt.plot(indices[model][0], indices[model][1])
            print(model)










