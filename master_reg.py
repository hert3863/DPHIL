#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 16:58:21 2018

@author: bakerh
"""

import numpy as np
from scipy import stats


def master_compare(exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    import glob
    import fnmatch
    both = {}

    for exp in exps:
        # import lists of successful files for each exp
        va = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                       exp + '/mon/va/*')
        tas = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                        exp + '/mon/tas/*')
        '''
        pr = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                       exp + '/mon/pr/*')
        ua = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                       exp + '/mon/ua/*')
        if exps[0] == 'Plus15-Future_LCO2':
            rlut = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                             exp + '/mon/rlut/*')
        '''
        # turn lists into just lists of exp IDs
        for i, item in enumerate(va):
            va[i] = va[i][79+len(exp):83+len(exp)]
        for i, item in enumerate(tas):
            tas[i] = tas[i][79+len(exp):83+len(exp)]
        '''
        for i, item in enumerate(ua):
            ua[i] = ua[i][79+len(exp):83+len(exp)]
        for i, item in enumerate(pr):
            pr[i] = pr[i][76+len(exp):80+len(exp)]
        if exps[0] == 'Plus15-Future_LCO2':
            for i, item in enumerate(rlut):
                rlut[i] = rlut[i][80+len(exp):84+len(exp)]
        '''
        b = []
        # compare lists and add to dictionary if both exist
        '''
        if exps[0] == 'Plus15-Future_LCO2':
            for i, item in enumerate(pr):
                if fnmatch.filter(tas, pr[i]) != [] and fnmatch.filter(ua, pr[i]) != [] and fnmatch.filter(va, pr[i]) != [] and fnmatch.filter(rlut, pr[i]) != []:
                    b.append(pr[i])
        else:
            for i, item in enumerate(pr):
                if fnmatch.filter(tas, pr[i]) != [] and fnmatch.filter(ua, pr[i]) != [] and fnmatch.filter(va, pr[i]) != []:
                    b.append(pr[i])
        '''
        for i, item in enumerate(tas):
            if fnmatch.filter(va, tas[i]) != []:
                b.append(tas[i])
        both[exp] = b
    return both


def master_load_had(both, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    from netCDF4 import Dataset

    def shiftlat(grid):
        latold = np.arange(89.375, -90.625, -1.25)
        latnew = np.arange(90, -91.25, -1.25)
        regrid = np.zeros((np.ma.size(grid, axis=0), 145, 192))
        for i in range(143):
            regrid[:, i+1, :] = ((grid[:, i, :]*np.cos(latold[i]*np.pi/180) +
                                  grid[:, i+1, :]*np.cos(latold[i+1]*np.pi/180)) /
                                 (2*np.cos(latnew[i+1]*np.pi/180)))
        return regrid

    va = {}
    ua = {}
    tas = {}
    pr = {}
    rlut = {}
    if exps[0] == 'Plus15-Future_LCO2':
        finstr = '_2090-01_2100-12'
    else:
        finstr = '_2005-01_2016-12'

    for exp in exps:
        a = []
        b = []
        c = []
        d = []
        e = []
        for i in range(np.ma.size(both[exp])):
            a.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/va/item15202_monthly_mean_' +
                     both[exp][i] + finstr + '.nc')
            b.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/ua/item15201_monthly_mean_' +
                     both[exp][i] + finstr + '.nc')
            c.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/tas/item3236_monthly_mean_' +
                     both[exp][i] + finstr + '.nc')
            d.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/pr/item5216_daily_mean_' +
                     both[exp][i] + finstr + '.nc.monmean.nc')
            if exps[0] == 'Plus15-Future_LCO2':
                e.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                         '/mon/rlut/item2205_monthly_mean_' +
                         both[exp][i] + finstr + '.nc')

        vabar = np.zeros((np.ma.size(a), 132, 145, 192))
        uabar = np.zeros((np.ma.size(a), 132, 145, 192))
        tasbar = np.zeros((np.ma.size(a), 132, 145, 192))
        prbar = np.zeros((np.ma.size(a), 132, 145, 192))
        rlutbar = np.zeros((np.ma.size(a), 132, 145, 192))

        for k, l in enumerate(a):
            nc_fid = Dataset(l, 'r')
            vabar[k] = shiftlat(nc_fid.variables['item15202_monthly_mean'][:, 2])
            print('Done va ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(b):
            nc_fid = Dataset(l, 'r')
            uabar[k] = shiftlat(nc_fid.variables['item15201_monthly_mean'][:, 2])
            print('Done ua ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(c):
            nc_fid = Dataset(l, 'r')
            tasbar[k] = nc_fid.variables['item3236_monthly_mean'][:, 0]
            print('Done tas ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(d):
            nc_fid = Dataset(l, 'r')
            prbar[k] = nc_fid.variables['item5216_daily_mean'][:, 0]
            print('Done pr ' + str(exp) + ' ' + str(k+1))
        if exps[0] == 'Plus15-Future_LCO2':
            for k, l in enumerate(e):
                nc_fid = Dataset(l, 'r')
                rlutbar[k] = nc_fid.variables['item2205_monthly_mean'][:, 0]
                print('Done rlut ' + str(exp) + ' ' + str(k+1))
            rlutbar = np.transpose(np.reshape(rlutbar, (-1, 12, 145, 192)),
                                   (1, 0, 2, 3))
            rlut[exp] = rlutbar
        vabar = np.transpose(np.reshape(vabar, (-1, 12, 145, 192)),
                             (1, 0, 2, 3))
        uabar = np.transpose(np.reshape(uabar, (-1, 12, 145, 192)),
                             (1, 0, 2, 3))
        tasbar = np.transpose(np.reshape(tasbar, (-1, 12, 145, 192)),
                              (1, 0, 2, 3))
        prbar = np.transpose(np.reshape(prbar, (-1, 12, 145, 192)),
                             (1, 0, 2, 3))
        va[exp] = vabar
        ua[exp] = uabar
        tas[exp] = tasbar
        pr[exp] = prbar * 86400

    return va, ua, tas, pr, rlut


def load_had(both, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    from netCDF4 import Dataset

    def shiftlat(grid):
        latold = np.arange(89.375, -90.625, -1.25)
        latnew = np.arange(90, -91.25, -1.25)
        regrid = np.zeros((np.ma.size(grid, axis=0), 145, 192))
        for i in range(143):
            regrid[:, i+1, :] = ((grid[:, i, :]*np.cos(latold[i]*np.pi/180) +
                                  grid[:, i+1, :]*np.cos(latold[i+1]*np.pi/180)) /
                                 (2*np.cos(latnew[i+1]*np.pi/180)))
        return regrid

    va = {}
    tas = {}
    if exps[0] == 'Plus15-Future_LCO2':
        finstr = '_2090-01_2100-12'
    else:
        finstr = '_2006-01_2016-12'
    sindices = np.zeros((10))
    for i in range(10):
        sindices[i] = i*12
    sindices = sindices.astype(int)
    for exp in exps:
        a = []
        c = []
        for i in range(np.ma.size(both[exp])):
            a.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/va/item15202_monthly_mean_' +
                     both[exp][i] + finstr + '.nc')
            c.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/tas/item3236_monthly_mean_' +
                     both[exp][i] + finstr + '.nc')

        vabar = np.zeros((np.ma.size(a), 145, 192))
        tasbar = np.zeros((np.ma.size(a), 145, 192))

        for k, l in enumerate(a):
            nc_fid = Dataset(l, 'r')
            vabar[k] = np.mean(shiftlat(nc_fid.variables['item15202_monthly_mean'][sindices, 2]), axis=0)
            print('Done va ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(c):
            nc_fid = Dataset(l, 'r')
            tasbar[k] = np.mean(nc_fid.variables['item3236_monthly_mean'][sindices, 0], axis=0)
            print('Done tas ' + str(exp) + ' ' + str(k+1))

        va[exp] = vabar
        tas[exp] = tasbar

    return va, tas


def master_load_mir(exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    from netCDF4 import Dataset
    import glob
    va = {}
    ua = {}
    pr = {}

    for exp in exps:
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp +
                             '/mon/va/*'))
        b = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp +
                             '/mon/ua/*'))
        d = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp +
                             '/mon/pr/*'))

        vabar = np.zeros((np.ma.size(a), 120, 128, 256))
        uabar = np.zeros((np.ma.size(a), 120, 128, 256))
        prbar = np.zeros((np.ma.size(a), 120, 128, 256))

        for k, l in enumerate(a):
            nc_fid = Dataset(l, 'r')
            vabar[k] = nc_fid.variables['va'][:, 9]
            print('Done va ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(b):
            nc_fid = Dataset(l, 'r')
            uabar[k] = nc_fid.variables['ua'][:, 9]
            print('Done ua ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(d):
            nc_fid = Dataset(l, 'r')
            prbar[k] = nc_fid.variables['pr'][:]
            print('Done pr ' + str(exp) + ' ' + str(k+1))

        vabar = np.transpose(np.reshape(vabar, (-1, 12, 128, 256)),
                             (1, 0, 2, 3))
        uabar = np.transpose(np.reshape(uabar, (-1, 12, 128, 256)),
                             (1, 0, 2, 3))
        prbar = np.transpose(np.reshape(prbar, (-1, 12, 128, 256)),
                             (1, 0, 2, 3))
        va[exp] = vabar
        ua[exp] = uabar
        pr[exp] = prbar * 86400

    return va, ua, pr


def load_mir(exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    from netCDF4 import Dataset
    import glob
    va = {}
    tas = {}
    sindices = np.zeros((10))
    for i in range(10):
        sindices[i] = i*12
    sindices = sindices.astype(int)
    for exp in exps:
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp +
                             '/mon/va/*'))
        b = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/' + exp +
                             '/mon/tas/*'))

        vabar = np.zeros((np.ma.size(a), 128, 256))
        tasbar = np.zeros((np.ma.size(a), 128, 256))

        for k, l in enumerate(a):
            nc_fid = Dataset(l, 'r')
            vabar[k] = np.mean(nc_fid.variables['va'][sindices, 9], axis=0)
            print('Done va ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(b):
            nc_fid = Dataset(l, 'r')
            tasbar[k] = np.mean(nc_fid.variables['tas'][sindices], axis=0)
            print('Done tas ' + str(exp) + ' ' + str(k+1))

        va[exp] = vabar
        tas[exp] = tasbar

    return va, tas


def master_load_cam(exps=['Plus15-Future_LCO2', 'Plus15-Future']):
    from netCDF4 import Dataset
    import glob
    va = {}
    ua = {}
    pr = {}

    for exp in exps:
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                             '/mon/va/*'))
        b = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                             '/mon/ua/*'))
        d = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                             '/mon/pr/*'))

        vabar = np.zeros((np.ma.size(a), 120, 96, 144))
        uabar = np.zeros((np.ma.size(a), 120, 96, 144))
        prbar = np.zeros((np.ma.size(a), 120, 96, 144))

        for k, l in enumerate(a):
            nc_fid = Dataset(l, 'r')
            vabar[k] = nc_fid.variables['va'][:, 9]
            print('Done va ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(b):
            nc_fid = Dataset(l, 'r')
            uabar[k] = nc_fid.variables['ua'][:, 9]
            print('Done ua ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(d):
            nc_fid = Dataset(l, 'r')
            prbar[k] = nc_fid.variables['pr'][:, 0]
            print('Done pr ' + str(exp) + ' ' + str(k+1))

        vabar = np.transpose(np.reshape(vabar, (-1, 12, 96, 144)),
                             (1, 0, 2, 3))
        uabar = np.transpose(np.reshape(uabar, (-1, 12, 96, 144)),
                             (1, 0, 2, 3))
        prbar = np.transpose(np.reshape(prbar, (-1, 12, 96, 144)),
                             (1, 0, 2, 3))
        va[exp] = vabar
        ua[exp] = uabar
        pr[exp] = prbar * 86400

    return va, ua, pr


def load_cam(exps=['Plus15-Future_LCO2', 'Plus15-Future']):
    from netCDF4 import Dataset
    import glob
    va = {}
    tas = {}
    sindices = np.zeros((10))
    for i in range(10):
        sindices[i] = i*12
    sindices = sindices.astype(int)
    for exp in exps:
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                             '/mon/va/*'))
        b = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                             '/mon/tas/*'))

        vabar = np.zeros((np.ma.size(a), 96, 144))
        tasbar = np.zeros((np.ma.size(a), 96, 144))

        for k, l in enumerate(a):
            nc_fid = Dataset(l, 'r')
            vabar[k] = np.mean(nc_fid.variables['va'][sindices, 9], axis=0)
            print('Done va ' + str(exp) + ' ' + str(k+1))
        for k, l in enumerate(b):
            nc_fid = Dataset(l, 'r')
            tasbar[k] = np.mean(nc_fid.variables['tas'][sindices, 0], axis=0)
            print('Done tas ' + str(exp) + ' ' + str(k+1))

        va[exp] = vabar
        tas[exp] = tasbar

    return va, tas


def master_regression_had(v, var, month_v, month_var,
                          exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2'],
                          mx=1, pre='no', ti=''):
    if month_var == 'JJA':
        var_l = var[exps[0]]
        var_h = var[exps[1]]
    else:
        var_l = var[exps[0]][month_var]
        var_h = var[exps[1]][month_var]
    if month_v == 'JJA':
        v_l = v[exps[0]]
        v_h = v[exps[1]]
    else:
        v_l = v[exps[0]][month_v]
        v_h = v[exps[1]][month_v]

    v_anom = v_h - v_l.mean(axis=0)
    var_anom = var_h - var_l.mean(axis=0)

    if month_v == 5:
        circ_ind = (v_anom[:, 32, 2] - v_anom[:, 42, 22] + v_anom[:, 22, 53] -
                    v_anom[:, 34, 108] + v_anom[:, 39, 124] -
                    v_anom[:, 34, 141] + v_anom[:, 37, 158] -
                    v_anom[:, 31, 177])
    elif month_v == 6:
        circ_ind = (v_anom[:, 30, 187] - v_anom[:, 32, 21] +
                    v_anom[:, 31, 40] - v_anom[:, 37, 93] +
                    v_anom[:, 41, 117] - v_anom[:, 35, 136] +
                    v_anom[:, 34, 152] - v_anom[:, 31, 168])
    elif month_v == 7:
        circ_ind = (v_anom[:, 44, 180] - v_anom[:, 34, 23] +
                    v_anom[:, 32, 43] - v_anom[:, 36, 96] +
                    v_anom[:, 38, 122] - v_anom[:, 33, 138] +
                    v_anom[:, 33, 155] - v_anom[:, 34, 170])
    else:
        circ_ind = (v_anom[:, 30, 0] - v_anom[:, 32, 21] +
                    v_anom[:, 31, 42] - v_anom[:, 37, 94] +
                    v_anom[:, 39, 121] - v_anom[:, 34, 138] +
                    v_anom[:, 35, 154] - v_anom[:, 30, 172])

    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()

    beta = np.zeros((145, 192))
    for i in range(145):
        for j in range(192):
            beta[i, j] = stats.linregress(circ_ind[:], (var_anom[:, i, j]))[0]
    maplotr(beta, 'had', colormax=mx, precip=pre, title=ti)
    # return beta


def master_regression_mir(v, var, month_v, month_var,
                          exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2'],
                          mx=1, pre='no', ti=''):
    if month_var == 'JJA':
        var_l = var[exps[0]]
        var_h = var[exps[1]]
    else:
        var_l = var[exps[0]][month_var]
        var_h = var[exps[1]][month_var]
    if month_v == 'JJA':
        v_l = v[exps[0]]
        v_h = v[exps[1]]
    else:
        v_l = v[exps[0]][month_v]
        v_h = v[exps[1]][month_v]

    v_anom = v_h - v_l.mean(axis=0)
    var_anom = var_h - var_l.mean(axis=0)

    if month_v == 5:
        circ_ind = (v_anom[:, 102, 3] - v_anom[:, 103, 22] +
                    v_anom[:, 105, 51] - v_anom[:, 98, 144] +
                    v_anom[:, 95, 166] - v_anom[:, 98, 186] +
                    v_anom[:, 98, 210] - v_anom[:, 97, 235])
    elif month_v == 6:
        circ_ind = (v_anom[:, 101, 250] - v_anom[:, 97, 32] +
                    v_anom[:, 98, 56] - v_anom[:, 97, 146] +
                    v_anom[:, 97, 164] - v_anom[:, 98, 181] +
                    v_anom[:, 98, 203] - v_anom[:, 100, 220])
    else:
        circ_ind = (v_anom[:, 95, 252] - v_anom[:, 95, 34] +
                    v_anom[:, 106, 73] - v_anom[:, 97, 145] +
                    v_anom[:, 97, 165] - v_anom[:, 98, 185] +
                    v_anom[:, 99, 205] - v_anom[:, 102, 232])

    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()

    beta = np.zeros((128, 256))
    for i in range(128):
        for j in range(256):
            beta[i, j] = stats.linregress(circ_ind[:], (var_anom[:, i, j]))[0]
    maplotr(beta, 'mir', colormax=mx, precip=pre, title=ti)
    # return beta


def master_regression_cam(v, var, month_v, month_var,
                          exps=['Plus15-Future_LCO2', 'Plus15-Future'],
                          mx=1, pre='no', ti=''):
    if month_var == 'JJA':
        var_l = var[exps[0]]
        var_h = var[exps[1]]
    else:
        var_l = var[exps[0]][month_var]
        var_h = var[exps[1]][month_var]
    if month_v == 'JJA':
        v_l = v[exps[0]]
        v_h = v[exps[1]]
    else:
        v_l = v[exps[0]][month_v]
        v_h = v[exps[1]][month_v]

    v_anom = v_h - v_l.mean(axis=0)
    var_anom = var_h - var_l.mean(axis=0)

    if month_v == 5:
        circ_ind = (v_anom[:, 73, 143] - v_anom[:, 72, 17] +
                    v_anom[:, 72, 32] + v_anom[:, 71, 92] -
                    v_anom[:, 71, 105] +
                    v_anom[:, 71, 116] - v_anom[:, 71, 127])
    elif month_v == 6:
        circ_ind = (v_anom[:, 73, 143] - v_anom[:, 72, 17] +
                    v_anom[:, 72, 32] + v_anom[:, 71, 92] -
                    v_anom[:, 71, 105] +
                    v_anom[:, 71, 116] - v_anom[:, 71, 127])
    elif month_v == 7:
        circ_ind = (v_anom[:, 73, 143] - v_anom[:, 72, 17] +
                    v_anom[:, 72, 32] + v_anom[:, 71, 92] -
                    v_anom[:, 71, 105] +
                    v_anom[:, 71, 116] - v_anom[:, 71, 127])
    else:
        circ_ind = (v_anom[:, 73, 143] - v_anom[:, 72, 17] +
                    v_anom[:, 72, 32] + v_anom[:, 71, 92] -
                    v_anom[:, 71, 105] +
                    v_anom[:, 71, 116] - v_anom[:, 71, 127])

    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()

    beta = np.zeros((96, 144))
    for i in range(96):
        for j in range(144):
            beta[i, j] = stats.linregress(circ_ind[:], (var_anom[:, i, j]))[0]
    maplotr(beta, 'cam', colormax=mx, precip=pre, title=ti)
    # return beta


def maplotr(pdata, model, colormax=1, colormin=-999, title='',
            precip='no'):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    if colormin == -999:
        colormin = -colormax
    if model == 'had':
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    if model == 'mir':
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
        lon = np.arange(0, 360, 360/256)
    if model == 'cam':
        lat = np.linspace(-90, 90, 96)
        lon = np.linspace(0, 357.5, 144)
    plt.figure()
    pdata, lon = shiftgrid(180., pdata, lon, start=False)
    pdata, lon = addcyclic(pdata, lon)
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
    ctrs = np.linspace(colormin, colormax, 17)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    b = plt.colorbar(plot, orientation='horizontal',
                     aspect=50, shrink=0.75, spacing='proportional')
    b.set_label(label=r'pr (mm day$^{-1}$)')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()
