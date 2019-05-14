# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 14:12:02 2017

@author: bakerh
"""
import numpy as np


def rws(uwnd, vwnd):
    from windspharm.standard import VectorWind

    uwnd = np.transpose(uwnd, (1, 2, 0))
    vwnd = np.transpose(vwnd, (1, 2, 0))
    w = VectorWind(uwnd, vwnd)
    eta = w.absolutevorticity()
    div = w.divergence()
    uchi, vchi = w.irrotationalcomponent()
    etax, etay = w.gradient(eta)
    s1 = -eta * div
    s2 = - (uchi * etax + vchi * etay)
    s = s1 + s2
    s = np.mean(s, axis=2)
    s1 = np.mean(s1, axis=2)
    s2 = np.mean(s2, axis=2)
    return s, s1, s2


def trop():
    from netCDF4 import Dataset
    from netcdfread import ncread
    import glob

    def lapserate(t, z17):
        """
        Produces plot of lapse rate of T data
        """
        import numpy as np
        k = 287.058 / 0.718
        L = np.zeros((17, 145))
        for i in range(16):
            L[i, :] = k * 9.81 * ((t[i+1, :] - t[i, :]) * (z17[i+1] + z17[i])) / (287.058*(t[i+1, :] + t[i, :]) * (z17[i+1] - z17[i]))
        return L

    sindices = np.zeros((11*3))
    windices = np.zeros((11*3))
    for i in range(11):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(11):
        windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    z = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/ta/item16203_monthly_mean_a011_2006-01_2016-12.nc', 'z9')
    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                  'All-Hist/mon/ta/*')
    tbar = np.zeros((np.ma.size(z), 145))

    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        tbar += np.mean(nc_fid.variables['item16203_monthly_mean'][sindices, :],
                        axis=(0, 3))
        print('Done ' + str(e+1))
    tbar /= np.ma.size(a)
    L = lapserate(tbar, z)
    return tbar, L


def monthlymeans(item, code, level, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11))
    for i in range(11):
        sindices[i] = i*12
    sindices = sindices.astype(int)
    output = {}

    exps = ['All-Nat', 'SST-Nat', 'GHG-Nat']

    for x, exp in enumerate(exps):
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                             exp + '/mon/' + item + '/*'))
        sdata = np.zeros((np.ma.size(a), 132, region[1]-region[0], region[3]-region[2]))
        tbar = np.zeros((12, region[1]-region[0], region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[e] = nc_fid.variables[code][:, level, region[0]:region[1],
                                              region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for i in range(12):
            tbar[i] = np.mean(sdata[:, sindices+i], axis=(0, 1))
        output[exp] = tbar
    return output


def dailymeans(item, code, level, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*30))
    for i in range(11):
        sindices[30*i:30*(i+1)] = np.arange(0, 30, 1) + 360*i
    sindices = sindices.astype(int)
    output = {}

    exps = ['All-Nat', 'SST-Nat', 'GHG-Nat']

    for x, exp in enumerate(exps):
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                             exp + '/day/' + item + '/*'))
        sdata = np.zeros((np.ma.size(a), 3960, region[1]-region[0], region[3]-region[2]))
        tbar = np.zeros((12, region[1]-region[0], region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[e] = nc_fid.variables[code][:, level, region[0]:region[1],
                                              region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for i in range(12):
            tbar[i] = np.mean(sdata[:, sindices+30*i], axis=(0, 1))
        output[exp] = tbar * 86400
    return output


def main_circ(v, f, field='pr', exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2'], lev=2, wind='no', daily='no', dfmax=2, r=[0, 145, 0, 192]):
    #b = compare_var(field, code, daily, exps)
    #v, f = importmeans(field, code, b, daily, wind, lev, exps)
    v521 = v[exps[0]][1]
    v522 = v[exps[1]][1]
    f521 = f[exps[0]][1]
    f522 = f[exps[1]][1]
    v521 = shiftlat(v521)
    v522 = shiftlat(v522)
    if wind == 'yes':
        f521 = shiftlat(f521)
        f522 = shiftlat(f522)
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v521 - v521.mean(axis=0)
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f521 - f521.mean(axis=0)
    '''
    v_circ = circ_pattern(d_v, v_anom, v_anom, r)
    v_res = d_v - v_circ
    maplot(d_v, 2, title='d_v')
    maplot(v_circ, 2, title='v_circ')
    maplot(v_res, 2, title='v_res')
    '''
    f_circ = circ_pattern(d_v, v_anom, f_anom, r)
    f_res = d_f - f_circ
    if field == 'pr':
        maplot(d_f, dfmax, title='d_' + field, precip='yes')
        maplot(f_circ, dfmax, title=field + '_circ', precip='yes')
        maplot(f_res, dfmax, title=field + '_residual', precip='yes')
    else:
        maplot(d_f, dfmax, title='d_' + field)
        maplot(f_circ, dfmax, title=field + '_circ')
        maplot(f_res, dfmax, title=field + '_residual')


'''
def main_circ_project(field, wind='no', daily='no', r=[0, 145, 0, 192]):
    import fnmatch
    b = compare_circ(field)
    b1 = compare_circ('item15202_monthly_mean')
    both = []
    for i, item in enumerate(b1):
        if fnmatch.filter(b, b1[i]) != []:
            both.append(b1[i])
    bth = {}
    bth['batch_521'] = both
    bth['batch_522'] = both
    v, f = importmeans(field, bth, daily, wind)
    v521 = v['batch_521'][1]
    v522 = v['batch_522'][1]
    f521 = f['batch_521'][1]
    f522 = f['batch_522'][1]
    v521 = shiftlat(v521)
    v522 = shiftlat(v522)
    if wind == 'yes':
        f521 = shiftlat(f521)
        f522 = shiftlat(f522)
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v522 - v521
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f522 - f521

    def circ_corr(d_v, v_anom, r=[0, 145, 0, 192]):
        lat = np.arange(90, -91.25, -1.25)
        meshlat = np.zeros((192, 145))
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
    maplot(d_v, 2, title='d_v')
    maplot(v_circ, 2, title='v_circ')
    maplot(v_res, 2, title='v_res')

    f_circ = np.sum(np.moveaxis(f_anom, 0, 2)*v_coef, axis=2)
    f_res = d_f - f_circ
    if field == 'item5216_daily_mean':
        maplot(d_f, 2, title='d_f', precip='yes')
        maplot(f_circ, 2, title='f_circ', precip='yes')
        maplot(f_res, 2, title='f_res', precip='yes')
    else:
        maplot(d_f, 2, title='d_f')
        maplot(f_circ, 2, title='f_circ')
        maplot(f_res, 2, title='f_res')
'''


def main_circ_reg(field, code, mxx=2, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2'], wind='no', lev=2, daily='no'):
    from scipy import stats
    bth = compare_var(field, code, daily, exps)
    months = [4]
    for month in months:
        v, f = importmeans(field, code, bth, daily, wind, lev, exps, month)
        v521 = v[exps[0]][1]
        v522 = v[exps[1]][1]
        f521 = f[exps[0]][1]
        f522 = f[exps[1]][1]
        v521 = shiftlat(v521)
        v522 = shiftlat(v522)
        if wind == 'yes':
            f521 = shiftlat(f521)
            f522 = shiftlat(f522)
        d_v = v522.mean(axis=0) - v521.mean(axis=0)
        v_anom = v522 - v521.mean(axis=0)
        d_f = f522.mean(axis=0) - f521.mean(axis=0)
        f_anom = f522 - f521.mean(axis=0)
    
        # circ_ind0 = (d_v[31, 0] - d_v[33, 21] + d_v[31, 42] - d_v[36, 72] +
                  #   d_v[39, 122] - d_v[34, 138] + d_v[35, 154] - d_v[31, 171])
        circ_ind = (v_anom[:, 31, 0] - v_anom[:, 33, 21] + v_anom[:, 31, 42] + v_anom[:, 39, 122] - v_anom[:, 34, 138] +
                    v_anom[:, 35, 154] - v_anom[:, 31, 171])  # - circ_ind0
        circ_ind = -1*v_anom[:, 34, 138]
        # f_ind0 = d_f[36, 10] + d_f[26, 52] + d_f[30, 135]
        # f_ind = f_anom[:, 36, 10] + f_anom[:, 26, 52] + f_anom[:, 30, 135] - f_ind0
        circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()
        beta = np.zeros((145, 192))
        for i in range(145):
            for j in range(192):
                beta[i, j] = stats.linregress(circ_ind[:], f_anom[:, i, j])[0]
        f_circ = beta
        f_res = d_f - f_circ
        if field == 'pr':
            #maplot(d_f, mxx, title=str(month) + ' d_' + field, precip='yes')
            maplot(f_circ, mxx, title=str(month) + ' ' + field + '_circ', precip='yes')
            #maplot(f_res, mxx, title=str(month) + ' ' + field + '_residual', precip='yes')
        else:
            maplot(d_f, mxx, title='d_' + field)
            maplot(f_circ, mxx, title=field + '_circ')
            maplot(f_res, mxx, title=field + '_residual')
    #return d_f, f_circ


def reg_code(v, f, exps, wind='no', precip='no', mxx=1):
    from scipy import stats
    v521 = v[exps[0]][1]
    v522 = v[exps[1]][1]
    f521 = f[exps[0]][1]
    f522 = f[exps[1]][1]
    v521 = shiftlat(v521)
    v522 = shiftlat(v522)
    if wind == 'yes':
        f521 = shiftlat(f521)
        f522 = shiftlat(f522)
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v522 - v521.mean(axis=0)
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f522 - f521.mean(axis=0)
    circ_ind = (v_anom[:, 31, 0] - v_anom[:, 33, 21] + v_anom[:, 31, 42] + v_anom[:, 39, 122] - v_anom[:, 34, 138] +
                v_anom[:, 35, 154] - v_anom[:, 31, 171])
    #circ_ind = v_anom[:, 39, 122] - v_anom[:, 34, 138] + v_anom[:, 35, 154]
    circ_ind = (v_anom[:, 30, 190] - v_anom[:, 33, 20] + v_anom[:, 31, 40] + v_anom[:, 41, 117] - v_anom[:, 36, 137] +
                v_anom[:, 34, 152] - v_anom[:, 31, 168])
    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()
    beta = np.zeros((145, 192))
    for i in range(145):
        for j in range(192):
            beta[i, j] = stats.linregress(circ_ind[:], (f_anom[:, i, j]-np.mean(f_anom[:, i, j])))[0]
    f_circ = beta
    f_res = d_f - f_circ
    maplot(f_circ, mxx, precip=precip)
    return f_circ



def main_circ_reg_orog(orog_x, orog_y, mxx=1e-3, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    from scipy import stats
    from scipy import interpolate
    bth = compare_var('uas', 'item3225_monthly_mean', 'no', exps)
    v, f = importmeans('uas', 'item3225_monthly_mean', bth, 'no', 'yes', 0, exps)
    v, f1 = importmeans('vas', 'item3226_monthly_mean', bth, 'no', 'yes', 0, exps)
    v521 = v[exps[0]][1]
    v522 = v[exps[1]][1]
    f521 = f[exps[0]][1]
    f522 = f[exps[1]][1]
    f521_1 = f1[exps[0]][1]
    f522_1 = f1[exps[1]][1]
    v521_ncep = np.zeros((np.ma.size(v521, axis=0), 73, 144))
    v522_ncep = np.zeros((np.ma.size(v522, axis=0), 73, 144))
    f521_ncep = np.zeros((np.ma.size(v521, axis=0), 73, 144))
    f522_ncep = np.zeros((np.ma.size(v522, axis=0), 73, 144))
    f521_1_ncep = np.zeros((np.ma.size(v521, axis=0), 73, 144))
    f522_1_ncep = np.zeros((np.ma.size(v522, axis=0), 73, 144))

    for i in range(np.ma.size(v521, axis=0)):
        h = interpolate.interp2d(lon, lat[::-1], v521[i, :, :])
        v521_ncep[i, :, :] = h(lon_ncep, lat_ncep[::-1])
        g = interpolate.interp2d(lon, lat[::-1], f521[i, :, :])
        f521_ncep[i, :, :] = g(lon_ncep, lat_ncep[::-1])
        g = interpolate.interp2d(lon, lat[::-1], f521_1[i, :, :])
        f521_1_ncep[i, :, :] = g(lon_ncep, lat_ncep[::-1])
    for i in range(np.ma.size(v522, axis=0)):
        h = interpolate.interp2d(lon, lat[::-1], v522[i, :, :])
        v522_ncep[i, :, :] = h(lon_ncep, lat_ncep[::-1])
        g = interpolate.interp2d(lon, lat[::-1], f522[i, :, :])
        f522_ncep[i, :, :] = g(lon_ncep, lat_ncep[::-1])
        g = interpolate.interp2d(lon, lat[::-1], f522_1[i, :, :])
        f522_1_ncep[i, :, :] = g(lon_ncep, lat_ncep[::-1])


    d_v = v522_ncep.mean(axis=0) - v521_ncep.mean(axis=0)
    v_anom = v522_ncep - v521_ncep.mean(axis=0)
    g521 = f521_ncep * orog_x + f521_1_ncep * orog_y
    g522 = f522_ncep * orog_x + f522_1_ncep * orog_y
    d_g = g522.mean(axis=0) - g521.mean(axis=0)
    g_anom = g522 - g521.mean(axis=0)

    #circ_ind = (v_anom[:, 31, 0] - v_anom[:, 33, 21] + v_anom[:, 31, 42] + v_anom[:, 39, 122] - v_anom[:, 34, 138] +
    #            v_anom[:, 35, 154] - v_anom[:, 31, 171])
    circ_ind = (np.mean(g_anom[:, 42:49, 114:121], axis=(1,2)))
    circ_ind = g_anom[:, 22, 29] - g_anom[:, 25, 33]

    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()
    beta = np.zeros((73, 144))
    for i in range(73):
        for j in range(144):
            beta[i, j] = stats.linregress(circ_ind[:], v_anom[:, i, j])[0]
    g_circ = beta
    g_res = d_v - g_circ
    maplot(d_v, lat_ncep, lon_ncep, mxx, title='d_v')
    maplot(g_circ, lat_ncep, lon_ncep, mxx, title='beta_v')
    maplot(g_res, lat_ncep, lon_ncep, mxx, title='omega_residual')
    return d_v, g_circ


def main_circ_reg_rw(field, code, mxx=2, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2'], wind='no', lev=2, daily='no'):
    from scipy import stats
    bth = compare_var(field, code, daily, exps)
    v, f = importmeans(field, code, bth, daily, wind, lev, exps)
    v521 = v[exps[0]][1]
    v522 = v[exps[1]][1]
    f521 = rws(f[exps[0]][1], v[exps[0]][1])[2]
    f522 = rws(f[exps[1]][1], v[exps[1]][1])[2]
    v521 = shiftlat(v521)
    v522 = shiftlat(v522)
    if wind == 'yes':
        f521 = shiftlat(f521)
        f522 = shiftlat(f522)
    d_v = v522.mean(axis=0) - v521.mean(axis=0)
    v_anom = v522 - v521.mean(axis=0)
    d_f = f522.mean(axis=0) - f521.mean(axis=0)
    f_anom = f522 - f521.mean(axis=0)

    # circ_ind0 = (d_v[31, 0] - d_v[33, 21] + d_v[31, 42] - d_v[36, 72] +
              #   d_v[39, 122] - d_v[34, 138] + d_v[35, 154] - d_v[31, 171])
    circ_ind = (v_anom[:, 31, 0] - v_anom[:, 33, 21] + v_anom[:, 31, 42] + v_anom[:, 39, 122] - v_anom[:, 34, 138] +
                v_anom[:, 35, 154] - v_anom[:, 31, 171])  # - circ_ind0
    # f_ind0 = d_f[36, 10] + d_f[26, 52] + d_f[30, 135]
    # f_ind = f_anom[:, 36, 10] + f_anom[:, 26, 52] + f_anom[:, 30, 135] - f_ind0
    circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()
    beta = np.zeros((145, 192))
    for i in range(145):
        for j in range(192):
            beta[i, j] = stats.linregress(circ_ind[:], f_anom[:, i, j])[0]
    f_circ = beta
    f_res = d_f - f_circ
    if field == 'pr':
        maplot(d_f, mxx, title='d_' + field, precip='yes')
        maplot(f_circ, mxx, title=field + '_circ', precip='yes')
        maplot(f_res, mxx, title=field + '_residual', precip='yes')
    else:
        maplot(d_f, mxx, title='d_' + field)
        maplot(f_circ, mxx, title=field + '_circ')
        maplot(f_res, mxx, title=field + '_residual')
    return d_f, f_circ


def main_circ_regT(field, code, wind='no', daily='no'):
    from scipy import stats
    import fnmatch
    b = compare_circ(field, code, daily)
    b1 = compare_circ('tas', 'item3236_monthly_mean', 'no')
    both = []
    for i, item in enumerate(b1):
        if fnmatch.filter(b, b1[i]) != []:
            both.append(b1[i])
    bth = {}
    bth['Plus15-Future_LCO2'] = both
    bth['Plus15-Future_HCO2'] = both
    v, f = importmeans(field, code, bth, daily, wind)
    v521 = v['Plus15-Future_LCO2'][1]
    v522 = v['Plus15-Future_HCO2'][1]
    f521 = f['Plus15-Future_LCO2'][1]
    f522 = f['Plus15-Future_HCO2'][1]
    if wind == 'yes':
        f521 = shiftlat(f521)
        f522 = shiftlat(f522)
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
    beta = np.zeros((145, 192))
    for i in range(145):
        for j in range(192):
            beta[i, j] = stats.linregress(circ_ind[:], f_anom[:, i, j])[0]
    f_circ = beta
    f_res = d_f - f_circ
    if field == 'pr':
        maplot(d_f, 2, title='d_f', precip='yes')
        maplot(f_circ, 2, title='f_circ', precip='yes')
        maplot(f_res, 2, title='f_res', precip='yes')
    else:
        maplot(d_f, 2, title='d_f')
        maplot(f_circ, 2, title='f_circ')
        maplot(f_res, 2, title='f_res')


def compare_var(items, code, daily, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2']):
    '''
    Compares the control and perturbed ensembles
    and outputs a list of exp IDs that have completed
    for both ensembles, coupled with the patch number
    '''
    import glob
    import fnmatch
    both = {}
    if daily == 'yes':
        time = 'day'
    else:
        time = 'mon'
    for exp in exps:
        # import lists of successful files for each exp
        t521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp + '/mon/va/*')
        i521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp + '/' + time + '/' + items + '/*')

        # turn lists into just lists of exp IDs
        for i, item in enumerate(t521):
            t521[i] = t521[i][55+len(exp)+len('va')+len('item15202_monthly_mean'):55+len(exp)+len('va')+len('item15202_monthly_mean')+4]
        for i, item in enumerate(i521):
            i521[i] = i521[i][55+len(exp)+len(items)+len(code):55+len(exp)+len(items)+len(code)+4]

        b = []
        # compare lists and add to dictionary if both exist
        for i, item in enumerate(i521):
            if fnmatch.filter(t521, i521[i]) != []:
                b.append(i521[i])
        both[exp] = b
    return both


def compare_circ(items, code, daily, exps):
    '''
    Compares the control and perturbed ensembles
    and outputs a list of exp IDs that have completed
    for both ensembles, coupled with the patch number
    '''
    import glob
    import fnmatch
    if daily == 'yes':
        time = 'day'
    else:
        time = 'mon'
    # import lists of successful files for each exp
    i521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exps[0] + '/' + time + '/' + items + '/*')
    i522 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exps[0] + '/' + time + '/' + items + '/*')

    # turn lists into just lists of exp IDs
    for i, item in enumerate(i522):
        i522[i] = i522[i][55+len(exps[1])+len(items)+len(code):59+len(exps[1])+len(items)+len(code)]
    for i, item in enumerate(i521):
        i521[i] = i521[i][55+len(exps[0])+len(items)+len(code):59+len(exps[0])+len(items)+len(code)]

    both = []
    # compare lists and add to dictionary if both exist
    for i, item in enumerate(i522):
        if fnmatch.filter(i521, i522[i]) != []:
            both.append(i522[i])
    return both


def importmeans(item, code, both, daily='no', wind='no', lev=2, exps=['Plus15-Future_LCO2', 'Plus15-Future_HCO2'], month='JJA'):
    from netCDF4 import Dataset
    if daily == 'yes':
        time = 'day'
    else:
        time = 'mon'
    sindicesm = np.zeros((11*3))
    windicesm = np.zeros((11*3))
    for i in range(11):
        sindicesm[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(11):
        windicesm[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindicesm = sindicesm.astype(int)
    windicesm = windicesm.astype(int)

    if month == 'JJA':
        m = 3
        n = 90
        sindices = np.zeros((11*3))
        windices = np.zeros((11*3))
        for i in range(11):
            sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
        for i in range(11):
            windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
        sindices = sindices.astype(int)
        windices = windices.astype(int)

        sindicesd = np.zeros((11*90))
        windicesd = np.zeros((11*90))
        for i in range(11):
            sindicesd[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
        for i in range(11):
            windicesd[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                      60), np.linspace(330+360*i,
                                                      359+360*i, 30)))
        sindicesd = sindicesd.astype(int)
        windicesd = windicesd.astype(int)
    else:
        m = 1
        n = 30
        sindices = np.zeros((11))
        windices = np.zeros((11))
        for i in range(11):
            sindices[i] = month+12*i
        for i in range(11):
            windices[i] = month+12*i
        sindices = sindices.astype(int)
        windices = windices.astype(int)

        sindicesd = np.zeros((11*30))
        windicesd = np.zeros((11*30))
        for i in range(11):
            sindicesd[30*i:30*(i+1)] = np.linspace(30*month+360*i, 30*(month+1)-1+360*i, 30)
        for i in range(11):
            windicesd[30*i:30*(i+1)] = np.linspace(30*month+360*i, 30*(month+1)-1+360*i, 30)
        sindicesd = sindicesd.astype(int)
        windicesd = windicesd.astype(int)

    outputv = {}
    outputf = {}

    for exp in exps:
        a = []
        for i in range(np.ma.size(both[exp])):
            a.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                     '/mon/va/item15202_monthly_mean_' +
                     both[exp][i] + '_2090-01_2100-12.nc')
        vbar = np.zeros((2, np.ma.size(a), 144, 192))
        sdata = np.zeros((3*11*np.ma.size(a), 144, 192))
        wdata = np.zeros((3*11*np.ma.size(a), 144, 192))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[3*11*e:3*11*(e+1)] = nc_fid.variables['item15202_monthly_mean'
                                                    ][sindicesm, 2, :]
            wdata[3*11*e:3*11*(e+1)] = nc_fid.variables['item15202_monthly_mean'
                                                    ][windicesm, 2, :]
            print('Done ' + str('va') + ' ' +
                  str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            vbar[0, y, :] = np.mean(wdata[3*11*y:3*11*(y+1), :], axis=0)
            vbar[1, y, :] = np.mean(sdata[3*11*y:3*11*(y+1), :], axis=0)
        outputv[exp] = vbar

    if wind == 'yes':
        le = 144
    elif item == 'zg' and daily == 'no':
        le = 145
        lev = 1
    else:
        le = 145
        lev = 0
    if daily == 'yes':
        for exp in exps:
            a = []
            for i in range(np.ma.size(both[exp])):
                a.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                         '/' + time + '/' + item + '/' + code + '_' +
                         both[exp][i])
            rbar = np.zeros((2, np.ma.size(a), le, 192))
            sdatad = np.zeros((n*11*np.ma.size(a), le, 192))
            wdatad = np.zeros((n*11*np.ma.size(a), le, 192))

            for e, d in enumerate(a):
                nc_fid = Dataset(d, 'r')
                sdatad[n*11*e:n*11*(e+1)] = nc_fid.variables[code
                                                           ][sindicesd, lev, :]
                wdatad[n*11*e:n*11*(e+1)] = nc_fid.variables[code
                                                           ][windicesd, lev, :]
                print('Done ' + str(item) + ' ' +
                      str(exp) + ' ' + str(e+1))
            for y in range(np.ma.size(a)):
                rbar[0, y, :] = np.mean(wdatad[n*11*y:n*11*(y+1), :], axis=0)
                rbar[1, y, :] = np.mean(sdatad[n*11*y:n*11*(y+1), :], axis=0)

            outputf[exp] = rbar * 86400
    else:
        for exp in exps:
            a = []
            for i in range(np.ma.size(both[exp])):
                a.append('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                         '/' + time + '/' + item + '/' + code + '_' +
                         both[exp][i] + '_2090-01_2100-12.nc')
            fbar = np.zeros((2, np.ma.size(a), le, 192))
            sdata = np.zeros((m*11*np.ma.size(a), le, 192))
            wdata = np.zeros((m*11*np.ma.size(a), le, 192))
            for e, d in enumerate(a):
                nc_fid = Dataset(d, 'r')
                sdata[m*11*e:m*11*(e+1)] = nc_fid.variables[code
                                                        ][sindices, lev, :]
                wdata[m*11*e:m*11*(e+1)] = nc_fid.variables[code
                                                        ][windices, lev, :]
                print('Done ' + str(item) + ' ' +
                      str(exp) + ' ' + str(e+1))
            for y in range(np.ma.size(a)):
                fbar[0, y, :] = np.mean(wdata[m*11*y:m*11*(y+1), :], axis=0)
                fbar[1, y, :] = np.mean(sdata[m*11*y:m*11*(y+1), :], axis=0)
            outputf[exp] = fbar
    return outputv, outputf


def shiftlat(grid):
    latold = np.arange(89.375, -90.625, -1.25)
    latnew = np.arange(90, -91.25, -1.25)
    regrid = np.zeros((np.ma.size(grid, axis=0), 145, 192))
    for i in range(143):
        regrid[:, i+1, :] = ((grid[:, i, :]*np.cos(latold[i]*np.pi/180) +
                              grid[:, i+1, :]*np.cos(latold[i+1]*np.pi/180)) /
                             (2*np.cos(latnew[i+1]*np.pi/180)))
    return regrid


def circ_pattern(d_v, v_anom, f_anom, r=[0, 145, 0, 192]):
    from statsmodels.api import OLS
    lat = np.arange(90, -91.25, -1.25)
    meshlat = np.zeros((192, 145))
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
    field_circ = np.zeros((np.ma.size(v_anom, axis=0), 145, 192))

    for a in range(np.ma.size(v_anom, axis=0)):
        field_circ[a, :] = v_coef[a] * f_anom[a, :]
    field_circ = np.sum(field_circ, axis=0)
    return field_circ


def stat_ross(u, lat, lon):
    '''
    This treatment makes the regions of reflection (beta<0 and u_m>0)
    and regions of absorption (beta>0 and u_m<0) indistinguishable.
    These regions of absorption can be identified separately by
    plotting the easterlies.
    '''
    from netcdfread import ncread
    #lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #u = np.mean(u, axis=0)
    meshlon, meshlat = np.meshgrid(lon, lat)
    phi = meshlat*np.pi/180
    cs_phi = np.cos(phi)
    dphi = (lat[0]-lat[1])*np.pi/180
    omega = 7.29e-5
    a = 6.371e6
    u_m = u / (a*cs_phi)
    u_1opp = np.gradient(u_m*cs_phi**2, dphi, axis=0) / cs_phi
    u_2opp = np.gradient(u_1opp, dphi, axis=0) / cs_phi
    beta = (2*omega - u_2opp)*cs_phi**2/a
    ks = np.sqrt(a*beta/u_m)
    return ks


def sphr2cart(grid):
    from scipy import interpolate
    import numpy as np
    from netcdfread import ncread
    r = 6.371e6
    x_tot = 2*np.pi * r * np.cos(60*np.pi/180)
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lon_r = np.radians(lon)
    lat_r = np.radians(lat)
    meshlon, meshlat = np.meshgrid(lon_r, lat_r)
    meshx = r * meshlon / 2
    meshy = r * meshlat
    f = interpolate.interp2d(lon_r, lat_r, grid)
    x = meshx[0, :]
    y = meshy[:, 0]
    grid_cart = f(x*lon_r.max()/x.max(), y*lat_r.max()/y.max())
    return x, y, grid_cart


def conts(data, level):
    ys = np.zeros((192))
    for i in range(192):
        z = data[:73, i]
        y = np.where(np.diff(np.sign(z-level)))[0]
        y1 = y + 1
        ys[i] = y + (level-z[y])*(y1-y)/(z[y1]-z[y])
    return ys


def ys2ya(ys, y):
    ya = np.zeros((192))
    for i in range(192):
        ya[i] = y[int(np.floor(ys[i]))] + (ys[i]-np.floor(ys[i]))*(y[int(np.floor(ys[i]))+1]-y[int(np.floor(ys[i]))])
    return ya


def arc_length(x, y):
    npts = len(x)
    arc = np.sqrt((x[1] - x[0])**2 + (y[1] - y[0])**2)
    for k in range(1, npts):
        arc = arc + np.sqrt((x[k] - x[k-1])**2 + (y[k] - y[k-1])**2)
    return arc


def iso(grid, level):
    r = 6.371e6
    x, y, grid_cart = sphr2cart(grid)
    ys = conts(grid_cart, level)
    ya = ys2ya(ys, y)
    arc = arc_length(x, ya)
    m = arc / (2*np.pi * r * np.cos(60*np.pi/180))
    return m


def maplot(pdata, colormax=1, colormin=-999, mask='no', title='', precip='no'):
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
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')   
    plt.figure()
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        pdata = np.ma.masked_array(pdata, mask=np.logical_not(lsm))
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


def maplot_orog(pdata, u, v, lat, lon, title=''):
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

    plt.figure()
    pdata, lon1 = shiftgrid(180., pdata, lon, start=False)
    pdata, lon1 = addcyclic(pdata, lon1)
    u, lon1 = shiftgrid(180., u, lon, start=False)
    u, lon1 = addcyclic(u, lon1)
    v, lon = shiftgrid(180., v, lon, start=False)
    v, lon = addcyclic(v, lon)
    meshlon, meshlat = np.meshgrid(lon, lat)

    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    plot = m.contour(x, y, pdata,
                     np.array([-3,-2.5,-2,-1.5,-1,-.5,.5,1,1.5,2,2.5,3])*1e-2,
                     colors='k')
    #plt.quiver(x[::3, ::3], y[::3, ::3], u[::3, ::3], v[::3, ::3])
    b = plt.colorbar(plot, orientation='horizontal',
                     aspect=50, shrink=0.75, spacing='proportional')
    b.set_label(label=r'$\omega$ (Pa s$^{-1}$)')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def maplot_cgt(pdata, colormax=1, mask='no', title='', precip='no'):
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
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    plt.figure()
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        pdata = np.ma.masked_array(pdata, mask=np.logical_not(lsm))
    pdata, lon = shiftgrid(250., pdata, lon, start=False)
    pdata, lon = addcyclic(pdata, lon)
    meshlon, meshlat = np.meshgrid(lon, lat)

    m = Basemap(projection='cyl', llcrnrlat=0, urcrnrlat=70,
                llcrnrlon=-110, urcrnrlon=250, resolution='c')
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
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional', label='Difference (m)')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def maplot_c(pdata, pdata1, title='', precip='no'):
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
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    plt.figure()

    pdata, lon1 = shiftgrid(180., pdata, lon, start=False)
    pdata1, lon = shiftgrid(180., pdata1, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)

    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    plot = m.contour(x, y, pdata, [4, 5, 6], colors='k')
    plt.clabel(plot, inline=1, fontsize=10)
    plot1 = m.contour(x, y, pdata1, [4, 5, 6], colors='r', linestyles='--')
    plt.clabel(plot, inline=1, fontsize=10)
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def maplot_subs(va500, cmax=1, colormin=-999, precip='no'):
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
    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    def plotter(pdata, colormax=1, colormin=-999, title=''):
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
        #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')   
        if colormin == -999:
            colormin = -colormax
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

        plt.title(title, y=1)
        plt.show()
        return plot

    ax1 = fig.add_subplot(3, 3, 1)
    plotter(np.mean(va500['All-Hist'][1],axis=0)-273.15,colormax=cmax*40,title='All-Hist (multiply scale by 40)')
    ax2 = fig.add_subplot(3, 3, 2)
    plotter(np.mean(va500['All-Hist'][1], axis=0)-np.mean(va500['All-Nat'][1],axis=0),colormax=cmax,title='All-Hist - All-Nat')
    ax3 = fig.add_subplot(3, 3, 3)
    plotter(np.mean(va500['Plus15-Future'][1], axis=0)-np.mean(va500['All-Hist'][1],axis=0),colormax=cmax,title='Plus15-Future - All-Hist')
    ax4 = fig.add_subplot(3, 3, 4)
    plotter(np.mean(va500['All-Hist'][1], axis=0)-np.mean(va500['GHG-Nat'][1],axis=0),colormax=cmax,title='All-Hist - GHG-Nat')
    ax5 = fig.add_subplot(3, 3, 5)
    plotter(np.mean(va500['SST-Nat'][1], axis=0)-np.mean(va500['All-Nat'][1],axis=0),colormax=cmax,title='SST-Nat - All-Nat')
    ax6 = fig.add_subplot(3, 3, 6)
    plotter(np.mean(va500['Plus15-Future_HCO2'][1], axis=0)-np.mean(va500['Plus15-Future_LCO2'][1],axis=0),colormax=cmax,title='Plus15-Future_HCO2 - Plus15-Future_LCO2')
    ax7 = fig.add_subplot(3, 3, 7)
    plotter(np.mean(va500['All-Hist'][1], axis=0)-np.mean(va500['SST-Nat'][1],axis=0),colormax=cmax,title='All-Hist - SST-Nat')
    ax8 = fig.add_subplot(3, 3, 9)
    plotter(np.mean(va500['Plus15-Future_LCO2'][1], axis=0)-np.mean(va500['All-Hist'][1],axis=0),colormax=cmax,title='Plus15-Future_LCO2 - All-Hist')
    ax9 = fig.add_subplot(3, 3, 8)
    plot = plotter(np.mean(va500['GHG-Nat'][1], axis=0)-np.mean(va500['All-Nat'][1],axis=0),colormax=cmax,title='GHG-Nat - All-Nat')

    cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max')
    
    b.set_label(label='t200 difference ($^\circ$C)', size=20, fontsize=20, fontname='Arial')
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.95)


















































