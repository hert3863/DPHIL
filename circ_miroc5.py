# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 13:14:16 2017

@author: bakerh
"""


import numpy as np
from netcdfread import ncread


def stat_rossM(u):
    '''
    This treatment makes the regions of reflection (beta<0 and u_m>0)
    and regions of absorption (beta>0 and u_m<0) indistinguishable.
    These regions of absorption can be identified separately by
    plotting the easterlies.
    '''
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon = np.arange(0, 360, 360/256)
    meshlon, meshlat = np.meshgrid(lon, lat)
    phi = meshlat*np.pi/180
    cs_phi = np.cos(phi)
    dphi = 1.25*np.pi/180
    omega = 7.29e-5
    a = 6.371e6
    u_m = u / (a*cs_phi)
    u_1opp = np.gradient(u_m*cs_phi**2, dphi, axis=0) / cs_phi
    u_2opp = np.gradient(u_1opp, dphi, axis=0) / cs_phi
    beta = (2*omega - u_2opp)*cs_phi**2/a
    ks = np.sqrt(a*beta/u_m)
    return ks


def main_circM(field, wind='no', daily='no', dfmax=2, r=[0, 128, 0, 256]):
    v, f = importmeansM(field, daily, wind)
    v521 = v['plus15_lower'][1]
    v522 = v['plus15_higher'][1]
    f521 = f['plus15_lower'][1]
    f522 = f['plus15_higher'][1]

    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v521 - v521.mean(axis=0)
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f521 - f521.mean(axis=0)
    '''
    v_circ = circ_patternC(d_v, v_anom, v_anom, r)
    v_res = d_v - v_circ
    maplotM(d_v, 2, title='d_v')
    maplotM(v_circ, 2, title='v_circ')
    maplotM(v_res, 2, title='v_res')
    '''
    f_circ = circ_patternM(d_v, v_anom, f_anom, r)
    f_res = d_f - f_circ
    if field == 'pr':
        maplotM(d_f*86400, dfmax, title='d_f', precip='yes')
        maplotM(f_circ*86400, dfmax, title='f_circ', precip='yes')
        maplotM(f_res*86400, dfmax, title='f_res', precip='yes')
    else:
        maplotM(d_f, dfmax, title='d_f')
        maplotM(f_circ, dfmax, title='f_circ')
        maplotM(f_res, dfmax, title='f_res')


'''
def main_circ_project(field, wind='no', daily='no', r=[0, 128, 0, 256]):
    import fnmatch
    b = compare_circC(field)
    b1 = compare_circC('item15202_monthly_mean')
    both = []
    for i, item in enumerate(b1):
        if fnmatch.filter(b, b1[i]) != []:
            both.append(b1[i])
    bth = {}
    bth['plus15wrong'] = both
    bth['plus15'] = both
    v, f = importmeansC(field, bth, daily, wind)
    v521 = v['plus15wrong'][1]
    v522 = v['plus15'][1]
    f521 = f['plus15wrong'][1]
    f522 = f['plus15'][1]
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v522 - v521
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f522 - f521

    def circ_corr(d_v, v_anom, r=[0, 128, 0, 256]):
        lat = np.arange(90, -91.25, -1.25)
        meshlat = np.zeros((256, 128))
        meshlat[:, :] = lat
        meshlat = np.transpose(meshlat)
        mlw = (np.cos(meshlat * np.pi/180))

        v_coef = np.zeros((np.ma.size(v_anom, axis=0)))
        v_anom_w = v_anom * np.sqrt(mlw)
        dv = d_v * np.sqrt(mlw)
        dv_sub = dv[r[0]:r[1], r[2]:r[3]]
        v_anom_sub = v_anom_w[:, r[0]:r[1], r[2]:r[3]]
        for i in range(np.ma.size(v_anom, axis=0)):
            v_coef[i] = (np.sum(dv_sub * v_anom_sub[i, :]) /
                         np.sqrt(np.sum(dv_sub * dv_sub) *
                         np.sum(v_anom_sub[i, :] * v_anom_sub[i, :])))
        return v_coef
    v_coef = circ_corr(d_v, v_anom, r)
    v_coef = v_coef / np.sum(v_coef)
    v_circ = np.sum(np.moveaxis(v_anom, 0, 2)*v_coef, axis=2)
    v_res = d_v - v_circ
    maplotM(d_v, 2, title='d_v')
    maplotM(v_circ, 2, title='v_circ')
    maplotM(v_res, 2, title='v_res')

    f_circ = np.sum(np.moveaxis(f_anom, 0, 2)*v_coef, axis=2)
    f_res = d_f - f_circ
    if field == 'item5216_daily_mean':
        maplotM(d_f, 2, title='d_f', precip='yes')
        maplotM(f_circ, 2, title='f_circ', precip='yes')
        maplotM(f_res, 2, title='f_res', precip='yes')
    else:
        maplotM(d_f, 2, title='d_f')
        maplotM(f_circ, 2, title='f_circ')
        maplotM(f_res, 2, title='f_res')
'''


def main_circ_regM(field, mxx=1, wind='no', daily='no'):
    from scipy import stats
    v, f = importmeansM(field, daily, wind)
    v521 = v['Plus15-Future_LCO2'][1]
    v522 = v['Plus15-Future_HCO2'][1]
    f521 = f['Plus15-Future_LCO2'][1]
    f522 = f['Plus15-Future_HCO2'][1]
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v522 - v521
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f522 - f521

    circ_ind = (v_anom[:, 95, 252]-v_anom[:, 102, 232]+v_anom[:, 99, 205]-v_anom[:, 98, 185]+
                v_anom[:, 97, 165]-v_anom[:, 97, 145])
    circ_ind = (circ_ind / np.mean(circ_ind))
    beta = np.zeros((128, 256))
    for i in range(128):
        for j in range(256):
            beta[i, j] = stats.linregress(circ_ind[:], f_anom[:, i, j])[0]
    f_circ = beta
    f_res = d_f - f_circ
    if field == 'pr':
        maplotM(d_f, mxx, title='d_' + field, precip='yes')
        maploMt(f_circ, mxx, title=field + '_circ', precip='yes')
        maplotM(f_res, mxx, title=field + '_residual', precip='yes')
    else:
        maplotM(d_f, mxx, title='d_' + field)
        maplotM(f_circ, mxx, title=field + '_circ')
        maplotM(f_res, mxx, title=field + '_residual')


def main_circ_regTM(field, wind='no', daily='no'):
    from scipy import stats
    v, f = importmeansM(field, daily, wind)
    v521 = v['plus15_lower'][1]
    v522 = v['plus15_higher'][1]
    f521 = f['plus15_lower'][1]
    f522 = f['plus15_higher'][1]
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v522 - v521
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f522 - f521

    # circ_ind0 = (d_v[31, 0] - d_v[33, 21] + d_v[31, 42] - d_v[36, 72] +
              #   d_v[39, 122] - d_v[34, 138] + d_v[35, 154] - d_v[31, 171])
    circ_ind = (v_anom[:, 36, 10] - v_anom[:, 25, 53] + v_anom[:, 30, 135])  # - circ_ind0
    # f_ind0 = d_f[36, 10] + d_f[26, 52] + d_f[30, 135]
    # f_ind = f_anom[:, 36, 10] + f_anom[:, 26, 52] + f_anom[:, 30, 135] - f_ind0
    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()
    beta = np.zeros((128, 256))
    for i in range(128):
        for j in range(256):
            beta[i, j] = stats.linregress(circ_ind[:], f_anom[:, i, j])[0]
    f_circ = beta
    f_res = d_f - f_circ
    if field == 'pr':
        maplotM(d_f, 2, title='d_f', precip='yes')
        maplotM(f_circ, 2, title='f_circ', precip='yes')
        maplotM(f_res, 2, title='f_res', precip='yes')
    else:
        maplotM(d_f, 2, title='d_f')
        maplotM(f_circ, 2, title='f_circ')
        maplotM(f_res, 2, title='f_res')


def importmeansM(item, daily='no', wind='no'):
    from netCDF4 import Dataset
    import glob
    sindices = np.zeros((10*3))
    windices = np.zeros((10*3))
    for i in range(10):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(10):
        windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    sindicesd = np.zeros((10*92))
    windicesd = np.zeros((10*90))
    for i in range(10):
        sindicesd[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windicesd[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                  59), np.linspace(334+365*i,
                                                  364+365*i, 31)))
    sindicesd = sindicesd.astype(int)
    windicesd = windicesd.astype(int)

    outputv = {}
    outputf = {}
    exps =  ['Plus15-Future_LCO2', 'Plus15-Future_HCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp + '/mon/va/*')
        vbar = np.zeros((2, np.ma.size(a), 128, 256))
        sdata = np.zeros((3*10*np.ma.size(a), 128, 256))
        wdata = np.zeros((3*10*np.ma.size(a), 128, 256))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[30*e:30*(e+1)] = nc_fid.variables['va'
                                                    ][sindices, 9, :]
            wdata[30*e:30*(e+1)] = nc_fid.variables['va'
                                                    ][windices, 9, :]
            print('Done ' + str('va') + ' ' +
                  str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            vbar[0, y, :] = np.mean(wdata[30*y:30*(y+1), :], axis=0)
            vbar[1, y, :] = np.mean(sdata[30*y:30*(y+1), :], axis=0)
        outputv[exp] = vbar

    if wind == 'yes':
        lev = 9
    if daily == 'yes':
        for x, exp in enumerate(exps):
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp + '/day/' + item + '/*')
            rbar = np.zeros((2, np.ma.size(a), 128, 256))
            sdatad = np.zeros((92*10*np.ma.size(a), 128, 256))
            wdatad = np.zeros((90*10*np.ma.size(a), 128, 256))

            for e, d in enumerate(a):
                nc_fid = Dataset(d, 'r')
                sdatad[920*e:920*(e+1)] = nc_fid.variables[item
                                                           ][sindicesd, :]
                wdatad[900*e:900*(e+1)] = nc_fid.variables[item
                                                           ][windicesd, :]
                print('Done ' + str(item) + ' ' +
                      str(exp) + ' ' + str(e+1))
            for y in range(np.ma.size(a)):
                rbar[0, y, :] = np.mean(wdatad[900*y:900*(y+1), :], axis=0)
                rbar[1, y, :] = np.mean(sdatad[920*y:920*(y+1), :], axis=0)

            outputf[exp] = rbar
    else:
        for x, exp in enumerate(exps):
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp + '/mon/' + item + '/*')
            fbar = np.zeros((2, np.ma.size(a), 128, 256))
            sdata = np.zeros((3*10*np.ma.size(a), 128, 256))
            wdata = np.zeros((3*10*np.ma.size(a), 128, 256))
            for e, d in enumerate(a):
                nc_fid = Dataset(d, 'r')
                sdata[30*e:30*(e+1)] = nc_fid.variables[item
                                                        ][sindices, lev, :]
                wdata[30*e:30*(e+1)] = nc_fid.variables[item
                                                        ][windices, lev, :]
                print('Done ' + str(item) + ' ' +
                      str(exp) + ' ' + str(e+1))
            for y in range(np.ma.size(a)):
                fbar[0, y, :] = np.mean(wdata[30*y:30*(y+1), :], axis=0)
                fbar[1, y, :] = np.mean(sdata[30*y:30*(y+1), :], axis=0)
            outputf[exp] = fbar
    return outputv, outputf


def circ_patternM(d_v, v_anom, f_anom, r=[0, 128, 0, 256]):
    from statsmodels.api import OLS
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/present/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon = np.arange(0, 360, 360/256)
    meshlat = np.zeros((256, 128))
    meshlat[:, :] = lat
    meshlat = np.transpose(meshlat)
    mlw = (np.cos(meshlat * np.pi/180))

    v_coef = np.zeros((np.ma.size(v_anom, axis=0)))
    v_anom_w = v_anom * mlw
    dv = d_v * mlw
    dv_sub = np.ndarray.flatten(dv[r[0]:r[1], r[2]:r[3]])
    v_anom_sub = np.transpose(np.reshape(v_anom_w[:, r[0]:r[1], r[2]:r[3]],
                                         (np.ma.size(v_anom, axis=0),
                                          (r[1]-r[0])*(r[3]-r[2]))))
    model = OLS(dv_sub, v_anom_sub).fit()
    v_coef = model._results.params
    field_circ = np.zeros((np.ma.size(v_anom, axis=0), 128, 256))

    for a in range(np.ma.size(v_anom, axis=0)):
        field_circ[a, :] = v_coef[a] * f_anom[a, :]
    field_circ = np.sum(field_circ, axis=0)
    return field_circ


def maplotM(pdata, colormax=1, mask='no', title='', precip='no'):
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
    from netcdfread import ncread
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon = np.arange(0, 360, 360/256)
    plt.figure()
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/Python/lsm_n72.nc', 'lsm')
        pdata = np.ma.masked_array(pdata, mask=np.logical_not(lsm))
    pdata, lon = shiftgrid(180., pdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)

    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    if precip == 'yes':
        my_cmap = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-colormax, colormax, 17)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    b = plt.colorbar(plot, orientation='horizontal',
                     aspect=50, shrink=0.75, spacing='proportional')
    b.set_label(label=r'-$\mathbf{v}_\chi\cdot\nabla\zeta$ (10$^{-12}$ s$^{-2}$)')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()

























