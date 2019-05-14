# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 13:43:32 2017

@author: bakerh
"""
import numpy as np


def main():
    wbgt, wbgt95 = wbgt_calc()
    tmeans = monthmeanTall()
    tx = tmax()
    r95 = R95p()
    r = dailymean('pr')
    tbar = dailymean('tas')
    # total_r, total_t, globprcnts = totals()
    # prcnts, all_mean = percents
    tx90p = {}
    tx90p['batch_518'] = tx['batch_518'][2]
    tx90p['batch_520'] = tx['batch_520'][2]
    tx90p['batch_521'] = tx['batch_521'][2]
    tx90p['batch_522'] = tx['batch_522'][2]
    tmeans_tx90p = monthmeanTall(reg='yes', regitem='item3236_daily_maximum')
    tmeans_r95p = monthmeanTall(reg='yes', regitem='item5216_daily_mean')
    tx_reg = tmax(reg='yes')
    tx90p_reg = {}
    tx90p_reg['batch_518'] = tx_reg['batch_518'][2]
    tx90p_reg['batch_520'] = tx_reg['batch_520'][2]
    tx90p_reg['batch_521'] = tx_reg['batch_521'][2]
    tx90p_reg['batch_522'] = tx_reg['batch_522'][2]
    r95_reg = R95p(reg='yes')
    a_t, b_t, c_t = mlrcluster(tmeans_tx90p, tx90p_reg)
    a_w, b_w, c_w = mlrcluster(tmeans, wbgt95)
    a_r, b_r, c_r = mlrcluster(tmeans_r95p, r95_reg, msk='TROP')
    return tbar, tmeans, tx, r, r95, wbgt95, a_t, b_t, c_t, a_w, b_w, c_w, a_r, b_r, c_r


def compare(items, exp):
    '''
    Compares the control and perturbed ensembles
    and outputs a list of exp IDs that have completed
    for both ensembles, coupled with the patch number
    '''
    import glob
    import fnmatch
    # import lists of successful files for each exp
    t521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/atmos/item3236_monthly_mean/*')
    i521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/atmos/' + items + '/*')

    # turn lists into just lists of exp IDs
    for i, item in enumerate(t521):
        t521[i] = t521[i][96:100]
    for i, item in enumerate(i521):
        i521[i] = i521[i][98+2*(len(items)-22):102+2*(len(items)-22)]

    both = []
    # compare lists and add to dictionary if both exist
    for i, item in enumerate(t521):
        if fnmatch.filter(i521, t521[i]) != []:
            both.append(t521[i])
    return both


def importall():
    wbgt, wbgt95 = wbgt_calc()
    tmeans = monthmeanTall()
    tx = tmax()
    r95 = R95ptot()
    r = dailymean('item5216_daily_mean')
    tbar = monthmean('item3236_monthly_mean', 0)
    return tbar, tmeans, tx, r, r95, wbgt95


def monthmean(item, code, level, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*3))
    windices = np.zeros((11*3))
    for i in range(11):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(11):
        windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}

    exps = ['All-Nat', 'SST-Nat', 'GHG-Nat']#, 'Plus15-Future_HCO2', 'Plus15-Future_LCO2']
    exps = ['Plus15-Future_HCO2', 'Plus15-Future_LCO2']

    for x, exp in enumerate(exps):
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                             exp + '/mon/' + item + '/*'))
        tbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((3*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((3*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[33*e:33*(e+1)] = nc_fid.variables[code
                                                    ][sindices, level,
                                                      region[0]:
                                                      region[1],
                                                      region[2]:region[3]]
            wdata[33*e:33*(e+1)] = nc_fid.variables[code
                                                    ][windices, level,
                                                      region[0]:
                                                      region[1],
                                                      region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            tbar[0, y, :] = np.mean(wdata[33*y:33*(y+1), :], axis=0)
            tbar[1, y, :] = np.mean(sdata[33*y:33*(y+1), :], axis=0)
        output[exp] = tbar
    return output


def dailymean(item, code, level=0, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}

    exps = ['All-Hist', 'Plus15-Future']

    for x, exp in enumerate(exps):
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' + exp +
                      '/day/' + item + '/*'))
        rbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables[code
                                                      ][sindices, level, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[990*e:990*(e+1)] = nc_fid.variables[code
                                                      ][windices, level, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            rbar[0, y, :] = np.mean(wdata[990*y:990*(y+1), :], axis=0)
            rbar[1, y, :] = np.mean(sdata[990*y:990*(y+1), :], axis=0)

        output[exp] = rbar  # * 86400
    return output


def R95ptot(region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item5216_daily_mean/*')
    sdatar95p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    wdatar95p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatar95p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        wdatar95p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        print('Done R95p calc ' + str(e+1))
    r95 = np.ones((2, region[1]-region[0], region[3]-region[2]))
    for i in range(region[1]-region[0]):
        for j in range(region[3]-region[2]):
            wdatar95p_d = wdatar95p[:, i, j] * 86400
            sdatar95p_d = sdatar95p[:, i, j] * 86400
            wdatar95p_d1 = wdatar95p_d[wdatar95p_d > 1]
            sdatar95p_d1 = sdatar95p_d[sdatar95p_d > 1]
            if np.ma.size(wdatar95p_d1) != 0:
                r95[0, i, j] = np.percentile(wdatar95p_d1, 95)
            if np.ma.size(sdatar95p_d1) != 0:
                r95[1, i, j] = np.percentile(sdatar95p_d1, 95)

    outputr95p = {}
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item5216_daily_mean/*')

        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        r95p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                        region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            wdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            print('Done precip ' + str(exp) + ' ' + str(e+1))
        for m in range(np.ma.size(a)):
            r95p[0, m] = np.sum(wdata[990*m:990*(m+1), :] *
                                (wdata[990*m:990*(m+1), :] >
                                r95[0]).astype(int), axis=0) / np.sum(wdata[990*m:990*(m+1), :], axis=0)
            r95p[1, m] = np.sum(sdata[990*m:990*(m+1), :] *
                                (sdata[990*m:990*(m+1), :] >
                                r95[1]).astype(int), axis=0) / np.sum(sdata[990*m:990*(m+1), :], axis=0)
        outputr95p[exp] = r95p / 0.01  # annual precip change
    return outputr95p


def R95p(region=[0, 145, 0, 192], reg='no'):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item5216_daily_mean/*')
    sdatar95p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    wdatar95p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatar95p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        wdatar95p[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        print('Done R95p calc ' + str(e+1))
    r95 = np.ones((2, region[1]-region[0], region[3]-region[2]))
    for i in range(region[1]-region[0]):
        for j in range(region[3]-region[2]):
            wdatar95p_d = wdatar95p[:, i, j] * 86400
            sdatar95p_d = sdatar95p[:, i, j] * 86400
            wdatar95p_d1 = wdatar95p_d[wdatar95p_d > 1]
            sdatar95p_d1 = sdatar95p_d[sdatar95p_d > 1]
            if np.ma.size(wdatar95p_d1) != 0:
                r95[0, i, j] = np.percentile(wdatar95p_d1, 95)
            if np.ma.size(sdatar95p_d1) != 0:
                r95[1, i, j] = np.percentile(sdatar95p_d1, 95)

    outputr95p = {}
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        if reg == 'yes':
            a1 = compare('item5216_daily_mean', exp)
            a = []
            for i in range(len(a1)):
                if exp == 'batch_518':
                    a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/\
atmos/item5216_daily_mean/item5216_daily_mean_' + a1[i] + '_2006-01_2016-12.nc')
                else:
                    a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/\
atmos/item5216_daily_mean/item5216_daily_mean_' + a1[i] + '_2090-01_2100-12.nc')
        else:
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                          '/atmos/item5216_daily_mean/*')

        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        r95p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                        region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            wdata[990*e:990*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            print('Done precip ' + str(exp) + ' ' + str(e+1))

        truth = np.zeros((2, 90*11*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(11*np.ma.size(a)):
            truth[0, 90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                          r95[0, :]).astype(int)
            truth[1, 90*k:90*(k+1), :] = (sdata[90*k:90*(k+1), :] >
                                          r95[1, :]).astype(int)
        for m in range(np.ma.size(a)):
            r95p[0, m, :] = np.sum(truth[0, 990*m:990*(m+1), :], axis=0)
            r95p[1, m, :] = np.sum(truth[1, 990*m:990*(m+1), :], axis=0)
        outputr95p[exp] = r95p / (11)  # days per season
    return outputr95p


def tmax(region=[0, 145, 0, 192], reg='no'):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*90))
    windices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    for i in range(11):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i, 59+360*i,
                                                 60), np.linspace(330+360*i,
                                                 359+360*i, 30)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/batch_518' +
                  '/atmos/item3236_daily_maximum/*')
    sdatatx90p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    wdatatx90p = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatatx90p[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                       ][sindices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        wdatatx90p[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                       ][windices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        print('Done TX90p calc ' + str(e+1))
    days = np.arange(0, 990*np.ma.size(a), 90).astype(int)
    t90 = np.zeros((2, 90, region[1]-region[0],
                    region[3]-region[2]))
    for i in range(90):
        t90[0, i, :] = np.percentile(wdatatx90p[days+i, :], 90, axis=0)
        t90[1, i, :] = np.percentile(sdatatx90p[days+i, :], 90, axis=0)

    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    output = {}

    for x, exp in enumerate(exps):
        if reg == 'yes':
            a1 = compare('item3236_daily_maximum', exp)
            a = []
            for i in range(len(a1)):
                if exp == 'batch_518':
                    a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/\
atmos/item3236_daily_maximum/item3236_daily_maximum_' + a1[i] + '_2006-01_2016-12.nc')
                else:
                    a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/\
atmos/item3236_daily_maximum/item3236_daily_maximum_' + a1[i] + '_2090-01_2100-12.nc')
        else:
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                          '/atmos/item3236_daily_maximum/*')
        tx90p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                          region[3]-region[2]))
        txx = np.zeros((2, np.ma.size(a), 11, region[1]-region[0],
                        region[3]-region[2]))
        tx = np.zeros((2, np.ma.size(a), region[1]-region[0],
                       region[3]-region[2]))
        sdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*11*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[990*e:990*(e+1)] = nc_fid.variables['item3236_daily_maximum'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done txmax ' + str(exp) + ' ' + str(e+1))

        truth = np.zeros((2, 90*11*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(11*np.ma.size(a)):
            truth[0, 90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                          t90[0, :]).astype(int)
            truth[1, 90*k:90*(k+1), :] = (sdata[90*k:90*(k+1), :] >
                                          t90[1, :]).astype(int)
        for m in range(np.ma.size(a)):
            tx90p[0, m, :] = np.sum(truth[0, 990*m:990*(m+1), :], axis=0)
            tx90p[1, m, :] = np.sum(truth[1, 990*m:990*(m+1), :], axis=0)
        tx90p = tx90p / (11)  # days per season

        for y in range(np.ma.size(a)):
            tx[0, y, :] = np.mean(wdata[990*y:990*(y+1), :], axis=0)
            tx[1, y, :] = np.mean(sdata[990*y:990*(y+1), :], axis=0)
            for t in range(11):
                txx[0, y, t, :] = np.max(wdata[990*y+90*t:990*y+90*(t+1), :],
                                         axis=0)
                txx[1, y, t, :] = np.max(sdata[990*y+90*t:990*y+90*(t+1), :],
                                         axis=0)
        txx = np.mean(txx, axis=2)
        output[exp] = tx, txx, tx90p
    return output


def wbgt_calc():
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((11*3))
    windices = np.zeros((11*3))
    for i in range(11):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    for i in range(11):
        windices[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}
    output_95p = {}

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + 'batch_518' +
                  '/atmos/item3236_monthly_mean/*')
    wbgt = np.zeros((2, np.ma.size(a), 145, 192))
    wbgt_95p = np.zeros((2, np.ma.size(a), 145, 192))
    sdataT = np.zeros((3*11*np.ma.size(a), 145, 192))
    wdataT = np.zeros((3*11*np.ma.size(a), 145, 192))
    sdataH = np.zeros((3*11*np.ma.size(a), 145, 192))
    wdataH = np.zeros((3*11*np.ma.size(a), 145, 192))
    wbgts = np.zeros((3*11*np.ma.size(a), 145, 192))
    wbgtw = np.zeros((3*11*np.ma.size(a), 145, 192))
    for e, d in enumerate(a):
        b = d[:58]+str(45)+d[60:80]+str(45)+d[82:]
        nc_fid = Dataset(d, 'r')
        nc_fib = Dataset(b, 'r')

        sdataT[33*e:33*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                 ][sindices, 0, :]
        wdataT[33*e:33*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                 ][windices, 0, :]
        sdataH[33*e:33*(e+1)] = nc_fib.variables['item3245_monthly_mean'
                                                 ][sindices, 0, :]
        wdataH[33*e:33*(e+1)] = nc_fib.variables['item3245_monthly_mean'
                                                 ][windices, 0, :]
        print('Done batch_518 ' + str(e+1))

    Es = 6.112*np.exp(17.62*(sdataT-273.15)/(243.12+sdataT-273.15))
    Ew = 6.112*np.exp(17.62*(wdataT-273.15)/(243.12+wdataT-273.15))
    es = sdataH*Es/100
    ew = wdataH*Ew/100
    wbgts = 0.567*(sdataT-273.15)+0.393*es+3.94
    wbgtw = 0.567*(wdataT-273.15)+0.393*ew+3.94
    for y in range(np.ma.size(a)):
        wbgt[0, y, :] = np.mean(wbgtw[33*y:33*(y+1), :], axis=0)
        wbgt[1, y, :] = np.mean(wbgts[33*y:33*(y+1), :], axis=0)

    days = np.arange(0, 33*np.ma.size(a), 3).astype(int)
    months = np.arange(0, 33, 3).astype(int)
    wbgt95 = np.zeros((2, 33, 145, 192))
    for i in range(3):
        wbgt95[0, months+i, :] = np.percentile(wbgtw[days+i, :], 95, axis=0)
        wbgt95[1, months+i, :] = np.percentile(wbgts[days+i, :], 95, axis=0)

    for m in range(np.ma.size(a)):
        wbgt_95p[0, m] = np.sum((wbgtw[33*m:33*(m+1), :] - wbgt95[0]) *
                                (wbgtw[33*m:33*(m+1), :] >
                                wbgt95[0]).astype(int), axis=0)
        wbgt_95p[1, m] = np.sum((wbgts[33*m:33*(m+1), :] - wbgt95[1]) *
                                (wbgts[33*m:33*(m+1), :] >
                                wbgt95[1]).astype(int), axis=0)

    output['batch_518'] = wbgt
    output_95p['batch_518'] = wbgt_95p/11

    exps = ['batch_520', 'batch_521', 'batch_522']

    for x, ep in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + ep +
                      '/atmos/item3236_monthly_mean/*')
        wbgt = np.zeros((2, np.ma.size(a), 145, 192))
        wbgt_95p = np.zeros((2, np.ma.size(a), 145, 192))
        sdataT = np.zeros((3*11*np.ma.size(a), 145, 192))
        wdataT = np.zeros((3*11*np.ma.size(a), 145, 192))
        sdataH = np.zeros((3*11*np.ma.size(a), 145, 192))
        wdataH = np.zeros((3*11*np.ma.size(a), 145, 192))
        wbgts = np.zeros((3*11*np.ma.size(a), 145, 192))
        wbgtw = np.zeros((3*11*np.ma.size(a), 145, 192))
        for e, d in enumerate(a):
            b = d[:58]+str(45)+d[60:80]+str(45)+d[82:]
            nc_fid = Dataset(d, 'r')
            nc_fib = Dataset(b, 'r')

            sdataT[33*e:33*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                     ][sindices, 0, :]
            wdataT[33*e:33*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                     ][windices, 0, :]
            sdataH[33*e:33*(e+1)] = nc_fib.variables['item3245_monthly_mean'
                                                     ][sindices, 0, :]
            wdataH[33*e:33*(e+1)] = nc_fib.variables['item3245_monthly_mean'
                                                     ][windices, 0, :]
            print('Done ' + str(ep) + ' ' + str(e+1))

        Es = 6.112*np.exp(17.62*(sdataT-273.15)/(243.12+sdataT-273.15))
        Ew = 6.112*np.exp(17.62*(wdataT-273.15)/(243.12+wdataT-273.15))
        es = sdataH*Es/100
        ew = wdataH*Ew/100
        wbgts = 0.567*(sdataT-273.15)+0.393*es+3.94
        wbgtw = 0.567*(wdataT-273.15)+0.393*ew+3.94
        for y in range(np.ma.size(a)):
            wbgt[0, y, :] = np.mean(wbgtw[33*y:33*(y+1), :], axis=0)
            wbgt[1, y, :] = np.mean(wbgts[33*y:33*(y+1), :], axis=0)

        for m in range(np.ma.size(a)):
            wbgt_95p[0, m] = np.sum((wbgtw[33*m:33*(m+1), :] - wbgt95[0]) *
                                    (wbgtw[33*m:33*(m+1), :] >
                                    wbgt95[0]).astype(int), axis=0)
            wbgt_95p[1, m] = np.sum((wbgts[33*m:33*(m+1), :] - wbgt95[1]) *
                                    (wbgts[33*m:33*(m+1), :] >
                                    wbgt95[1]).astype(int), axis=0)
        output[ep] = wbgt
        output_95p[ep] = wbgt_95p/11
    return output, output_95p


def wbgt_test(t, td):
    def wbgt(t, td):
        rh = 100 * (np.exp((17.625*td)/(243.04+td)) /
                    np.exp((17.625*t)/(243.04+t)))
        E = 6.112*np.exp(17.62*(t)/(243.12+t))
        e = rh*E/100
        wbgt = 0.567*(t)+0.393*e+3.94
        return wbgt
    wbgt_sim = wbgt(t, td)
    wbgt_sim = np.mean(wbgt_sim)
    t_d = np.mean(np.reshape(t, (-1, 4)), axis=1)
    td_d = np.mean(np.reshape(td, (-1, 4)), axis=1)
    wbgt_dailymean = np.mean(wbgt(t_d, td_d))
    return wbgt_sim, wbgt_dailymean


def latweightmean(data, msk='no'):
    from netcdfread import ncread
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    weighted = data * meshlatweight
    if msk == 'yes':
        mask = np.zeros((np.ma.size(data, axis=0), 145, 192))
        mask[:] = lsm1
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'NH':
        mask = np.zeros((np.ma.size(data, axis=0), 73, 192))
        lsm1 = lsm1[:73, :]
        mask[:] = lsm1
        weighted = weighted[:, :73, :]
        meshlatweight = meshlatweight[:73, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'TROP':
        mask = np.zeros((np.ma.size(data, axis=0), 49, 192))
        lsm1 = lsm1[48:97, :]
        mask[:] = lsm1
        weighted = weighted[:, 48:97, :]
        meshlatweight = meshlatweight[48:97, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'MID':
        mask = np.zeros((np.ma.size(data, axis=0), 48, 192))
        lsm1 = lsm1[:48, :]
        mask[:] = lsm1
        weighted = weighted[:,  :48, :]
        meshlatweight = meshlatweight[:48, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'no':
        meaned = np.mean(weighted) / np.mean(meshlatweight)
    elif np.ma.size(msk) == 5:
        mask = np.zeros((np.ma.size(data, axis=0), msk[2]-msk[1],
                         192+msk[4]-msk[3]))
        lsm1 = np.concatenate((lsm1[msk[1]:msk[2], msk[3]:], lsm1[msk[1]:msk[2], :msk[4]]), axis=1)
        mask[:] = lsm1
        weighted = np.concatenate((weighted[:, msk[1]:msk[2], msk[3]:], weighted[:, msk[1]:msk[2], :msk[4]]), axis=2)
        meshlatweight = np.concatenate((meshlatweight[msk[1]:msk[2], msk[3]:], meshlatweight[msk[1]:msk[2], :msk[4]]), axis=1)
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    else:
        mask = np.zeros((np.ma.size(data, axis=0), msk[1]-msk[0],
                         msk[3]-msk[2]))
        lsm1 = lsm1[msk[0]:msk[1], msk[2]:msk[3]]
        mask[:] = lsm1
        weighted = weighted[:, msk[0]:msk[1], msk[2]:msk[3]]
        meshlatweight = meshlatweight[msk[0]:msk[1], msk[2]:msk[3]]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    return meaned


def totals():
    import glob
    from netCDF4 import Dataset
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    total_r = np.zeros((4, 2, 3))  # exp; global, land only;
    total_t = np.zeros((4, 2, 3))  # all, winter, summer
    percentiles = np.zeros((4, 4, 5))
    globprcnts = np.zeros((4, 4, 5))
    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item5216_daily_mean/*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                      '/atmos/item3236_monthly_mean/*')
        dataR = np.zeros((360*11*np.ma.size(a), 145, 192))
        dataT = np.zeros((12*11*np.ma.size(b), 145, 192))

        for e, d in enumerate(a):
            sindices = np.zeros((11*90*np.ma.size(a)))
            windices = np.zeros((11*90*np.ma.size(a)))
            for i in range(11*np.ma.size(a)):
                sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
            for i in range(11*np.ma.size(a)):
                windices[90*i:90*(i+1)] = np.concatenate((np.linspace(360*i,
                                                          59+360*i, 60),
                                                         np.linspace(330+360*i,
                                                         359+360*i, 30)))
            sindices = sindices.astype(int)
            windices = windices.astype(int)

            nc_fid = Dataset(d, 'r')
            dataR[3960*e:3960*(e+1)] = nc_fid.variables['item5216_daily_mean'
                                                        ][:, 0, :]
            print('Done R totals ' + str(exp) + ' ' + str(e+1))
        total_r[x, 0, 0] = latweightmean(dataR, msk='no') * 86400
        total_r[x, 0, 1] = latweightmean(dataR[windices, :], msk='no') * 86400
        total_r[x, 0, 2] = latweightmean(dataR[sindices, :], msk='no') * 86400
        total_r[x, 1, 0] = latweightmean(dataR, msk='yes') * 86400
        dataRw = dataR[windices, :]
        dataRs = dataR[sindices, :]
        total_r[x, 1, 1] = latweightmean(dataRw, msk='yes') * 86400
        total_r[x, 1, 2] = latweightmean(dataRs, msk='yes') * 86400
        '''
        dataR_point = np.zeros((np.ma.size(a)*3960))
        for l in range(np.ma.size(a)*3960):
            dataR_point[l] = latweightmean(dataR[l, :]) * 86400
        dataR_11 = np.zeros((np.ma.size(a)))
        for y in range(np.ma.size(a)):
            dataR_11[y] = np.mean(dataR_point[3960*y:3960*(y+1)], axis=0)
        '''
        dataR_11 = np.zeros((2, np.ma.size(a)))
        for y in range(np.ma.size(a)):
            dataR_11[0, y] = latweightmean(dataR[3960*y:3960*(y+1)], msk='no') * 86400
            dataR_11[1, y] = latweightmean(dataR[3960*y:3960*(y+1)], msk='yes') * 86400
        percentiles[2, x, :] = np.percentile(dataR_11[0, :], [10, 50, 90, 5, 95])
        percentiles[3, x, :] = np.percentile(dataR_11[1, :], [10, 50, 90, 5, 95])

        for f, g in enumerate(b):
            sindicesm = np.zeros((11*3*np.ma.size(b)))
            windicesm = np.zeros((11*3*np.ma.size(b)))
            for i in range(11*np.ma.size(b)):
                sindicesm[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
            for i in range(11*np.ma.size(b)):
                windicesm[3*i:3*(i+1)] = [12*i, 1+12*i, 11+12*i]
            sindicesm = sindicesm.astype(int)
            windicesm = windicesm.astype(int)
            nc_fid = Dataset(g, 'r')
            dataT[132*f:132*(f+1)] = nc_fid.variables['item3236_monthly_mean'
                                                      ][:, 0, :]
            print('Done T totals ' + str(exp) + ' ' + str(f+1))
        total_t[x, 0, 0] = latweightmean(dataT, msk='no')
        total_t[x, 0, 1] = latweightmean(dataT[windicesm, :], msk='no')
        total_t[x, 0, 2] = latweightmean(dataT[sindicesm, :], msk='no')
        total_t[x, 1, 0] = latweightmean(dataT, msk='yes')
        dataTw = dataT[windicesm, :]
        dataTs = dataT[sindicesm, :]
        total_t[x, 1, 1] = latweightmean(dataTw, msk='yes')
        total_t[x, 1, 2] = latweightmean(dataTs, msk='yes')
        '''
        dataT_point = np.zeros((np.ma.size(b)*132))
        for l in range(np.ma.size(b)*132):
            dataT_point[l] = latweightmean(dataT[l, :])
        dataT_11 = np.zeros((np.ma.size(b)))
        for y in range(np.ma.size(b)):
            dataT_11[y] = np.mean(dataT_point[132*y:132*(y+1)], axis=0)
        '''
        dataT_11 = np.zeros((2, np.ma.size(b)))
        for y in range(np.ma.size(b)):
            dataT_11[0, y] = latweightmean(dataT[132*y:132*(y+1)], msk='no')
            dataT_11[1, y] = latweightmean(dataT[132*y:132*(y+1)], msk='yes')
        percentiles[0, x, :] = np.percentile(dataT_11[0, :], [10, 50, 90, 5, 95])
        percentiles[1, x, :] = np.percentile(dataT_11[1, :], [10, 50, 90, 5, 95])
    for b in range(4):
        for a in range(4):
            globprcnts[b, a, 0] = percentiles[b, a, 1]
            globprcnts[b, a, 1] = (percentiles[b, a, 1] -
                                   percentiles[b, a, 0])
            globprcnts[b, a, 2] = (percentiles[b, a, 2] -
                                   percentiles[b, a, 1])
            globprcnts[b, a, 3] = (percentiles[b, a, 1] -
                                   percentiles[b, a, 3])
            globprcnts[b, a, 4] = (percentiles[b, a, 4] -
                                   percentiles[b, a, 1])
    return total_r, total_t, globprcnts


def monthmeanTall(reg='no', regitem='item5216_daily_mean'):
    import glob
    from netCDF4 import Dataset

    output = {}

    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        if reg == 'yes':
            a1 = compare(regitem, exp)
            a = []
            for i in range(len(a1)):
                if exp == 'batch_518':
                    a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/\
atmos/item3236_monthly_mean/item3236_monthly_mean_' + a1[i] + '_2006-01_2016-12.nc')
                else:
                    a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/\
atmos/item3236_monthly_mean/item3236_monthly_mean_' + a1[i] + '_2090-01_2100-12.nc')
        else:
            a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                          '/atmos/item3236_monthly_mean/*')
        data = np.zeros((132*np.ma.size(a), 145, 192))
        tbar = np.zeros((np.ma.size(a)))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            data[132*e:132*(e+1)] = nc_fid.variables['item3236_monthly_mean'
                                                     ][:, 0, :]
            print('Done item3236_monthly_mean ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            tbar[y] = latweightmean(data[132*y:132*(y+1), :], msk='no')
        output[exp] = tbar
    return output


def percents(tbar, txmax, wbgt_95p, r, r95p):
    exps = ['batch_518', 'batch_520', 'batch_521', 'batch_522']
    percentiles = np.zeros((4, 10, 3))
    prcnts = np.zeros((4, 10, 3))
    all_mean = np.zeros((4, 10))

    for e, exp in enumerate(exps):
        tbar1 = tbar[exp][1, :, :, :]
        tx90p1 = txmax[exp][2][1, :, :, :]
        wbgt_95p1 = wbgt_95p[exp][1, :, :, :]
        r1 = r[exp][1, :, :, :]
        r95p1 = r95p[exp][1, :, :, :]
        varis = [tbar1, tx90p1, wbgt_95p1, r1, r95p1]
        for i, var in enumerate(varis):
            means = []
            means_trop = []
            for j in range(np.ma.size(var, axis=0)):
                means.append(latweightmean(np.expand_dims(var[j, :],
                                                          axis=0), msk='MID'))
                means_trop.append(latweightmean(np.expand_dims(var[j, :],
                                                               axis=0),
                                                msk='TROP'))
            percentiles[e, 2*i, :] = np.percentile(means, [10, 50, 90])
            percentiles[e, 2*i+1, :] = np.percentile(means_trop, [10, 50, 90])
            all_mean[e, 2*i] = np.mean(means)
            all_mean[e, 2*i+1] = np.mean(means_trop)

    for a in range(4):
        for b in range(10):
            prcnts[a, b, 0] = percentiles[a, b, 1]
            prcnts[a, b, 1] = percentiles[a, b, 1] - percentiles[a, b, 0]
            prcnts[a, b, 2] = percentiles[a, b, 2] - percentiles[a, b, 1]
    return prcnts, all_mean


def mlr(tmean, var, norm='no'):
    import pandas as pd
    from statsmodels.formula.api import ols

    co2 = np.array([390.4, 423.4, 395.8, 550.0])
    f = 3.74 * np.log(co2/278) / np.log(2)  # 278 as preindust [CO2]
    d_f = f[1:] - f[0]
    d_t = tmean[1:]-tmean[0]
    d_var = var[1:]-var[0]
    # d_co2 = co2-co2[0]
    # logco2 = np.log(co2)
    # d_logco2 = logco2-logco2[0]
    if norm == 'yes':
        d_t /= d_t[1]
        d_var /= d_var[1]
    data = pd.DataFrame({'d_f': d_f, 'd_t': d_t, 'd_var': d_var})
    model = ols("d_var ~ d_t + d_f - 1", data).fit()
    # print(model.summary())
    # print(model._results.params)
    return model._results.params


def mlrcluster(tmeans, varmeans, msk='MID'):
    import pandas as pd
    from statsmodels.formula.api import ols

    def elevenyrmns(data, msk='MID'):
        exps = ['batch_520', 'batch_521', 'batch_522']
        pmean = []
        for j in range(np.ma.size(data['batch_518'][1], axis=0)):
            pmean.append(latweightmean(np.expand_dims(data['batch_518']
                                                          [1, j, :],
                                                      axis=0), msk))
        pmean = np.mean(pmean)
        means = []
        size = []
        for exp in exps:
            size.append(np.ma.size(data[exp][1], axis=0))
            for j in range(np.ma.size(data[exp][1], axis=0)):
                means.append(latweightmean(np.expand_dims(data[exp][1, j, :],
                                                          axis=0), msk)-pmean)
        return means, size
    co2 = np.array([390.4, 423.4, 395.8, 550.0])
    f = 3.74 * np.log(co2/278) / np.log(2)  # 278 as preindust [CO2]
    d_f = f[1:] - f[0]
    t_means = (np.concatenate((tmeans['batch_520'], tmeans['batch_521'],
                               tmeans['batch_522'])) -
               np.mean(tmeans['batch_518']))
    t_size = np.ma.size(t_means)
    var_means, v_size = elevenyrmns(varmeans, msk)
    if t_size == np.sum(v_size):
        d_fs = np.zeros((np.sum(v_size)))
        d_fs[:v_size[0]] = d_f[0]
        d_fs[v_size[0]:v_size[0] + v_size[1]] = d_f[1]
        d_fs[v_size[0] + v_size[1]:] = d_f[2]
    else:
        print("Error, Tmean and extreme differ in size")
    data = pd.DataFrame({'d_f': d_fs, 'd_t': t_means, 'd_var': var_means})
    model = ols("d_var ~ d_t + d_f - 1", data).fit()
    # print(model.summary())
    # print(model._results.params)
    beta, alpha = model._results.params
    cvs = model.cov_params()
    return alpha, beta, cvs


def fig4(tmeans, conversion, co2plume, rcpplume, alpha, beta, cvs, title=''):
    from matplotlib import pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    def gtc2ppm(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        return ["%.0f" % z for z in x2]

    def gtc2ppmfloat(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        return x2

    def ppm2gtc(x1):
        x2 = np.interp(x1, conversion[:, 1], conversion[:, 2])
        return x2

    def df2gtc(x1):
        dco2 = 390.4*np.exp(np.log(2)*x1/3.74)
        x2 = np.interp(dco2, conversion[:, 1], conversion[:, 2])
        return x2

    def gtc2df(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        x3 = (3.74/np.log(2))*np.log(x2/390.4)
        return x3

    def dfdc(c):
        ppm1 = gtc2ppmfloat(c-10)
        ppm2 = gtc2ppmfloat(c+10)
        f1 = 3.74 * np.log(ppm1/278) / np.log(2)
        f2 = 3.74 * np.log(ppm2/278) / np.log(2)
        dfdc = (f2-f1)/20
        return dfdc

    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    co2 = np.array([390.4, 423.4, 395.8, 550.0])
    tmean = np.zeros((4))
    tmean[0] = np.mean(tmeans['batch_518'])
    tmean[1] = np.mean(tmeans['batch_520'])
    tmean[2] = np.mean(tmeans['batch_521'])
    tmean[3] = np.mean(tmeans['batch_522'])

    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv['d_f'][1] +
                       (x*alpha/beta**2)**2*cv['d_t'][0] -
                       2*(x**2*alpha*cv['d_t'][1]**2/beta**3))
        return sigT

    co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301

    ax1.set_ylim([0, 1.5])
    hfont = {'fontname': 'Arial'}
    ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2016 \
($^{\circ}$C)', fontsize=20, **hfont)
    ax1.set_xlabel('Cumulative total anthropogenic CO$\mathregular{_{2}}$ \
emissions from 2011 (Gt C)', fontsize=20, **hfont)
    ax1.set_xticks(np.arange(0, 1100, 100))
    ax1.set_xlim([0, 650])

    x11 = rcpplume[0, 20]
    y11 = rcpplume[1, 20]
    x22 = rcpplume[0, 3]
    y22 = rcpplume[1, 3]
    g1 = (rcpplume[1, 19]-rcpplume[1, 20])/(rcpplume[0, 19]-rcpplume[0, 20])
    g2 = (rcpplume[1, 4]-rcpplume[1, 3])/(rcpplume[0, 4]-rcpplume[0, 3])
    ax1.axhline(tmean[1]-tmean[0], linestyle='--', color='k', linewidth=2,
                dashes=(10, 10))
    # ax1.fill_between(co2plume[0, :], co2plume[1, :], facecolor='gray',
    #                  alpha=0.2)
    ax1.fill_between(rcpplume[0, :], rcpplume[1, :], facecolor='#F5D9CE',
                     edgecolor='#eebdaa', alpha=1)
    x = np.arange(0, 1010, .01)
    y = (-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0])
    cns, = ax1.plot(x, y, linewidth=2, color='#0708EC')
    y1 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) -
          1.6445*sigt(gtc2df(x), cvs, alpha, beta))
    y2 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) +
          1.6445*sigt(gtc2df(x), cvs, alpha, beta))
    ax1.fill_between(x, y1, y2, facecolor='#0708EC', edgecolor='#0708EC',
                     alpha=.25)

    oub = (tmean[1]-tmean[0] - y22) / g2 + x22
    lb = (tmean[1]-tmean[0] - y11) / g1 + x11
    ub = x[np.argwhere(np.diff(np.sign(y -
                               g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubu = x[np.argwhere(np.diff(np.sign(y2 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubl = x[np.argwhere(np.diff(np.sign(y1 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    spr = ax1.axvspan(ubl, ubu, color='#ff6600', alpha=0.25)
    ax1.axvline(oub, linestyle='--', color='#FF0400', linewidth=2)
    ax1.axvline(lb, linestyle='--', color='#FF0400', linewidth=2)
    ubp = ax1.axvline(ub, linestyle='--', color='#FF7F2A', linewidth=2)
    mns = ax1.scatter(rcp_conc, tmean[1]-tmean[0], color='#5EB3FA',
                      marker='o', s=100, lw=3, label='RCP2.6 2090s')
    '''
    ax2 = ax1.twiny()
    ax2.spines['left'].set_position(('axes', 0.0))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 1.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax2.set_xlabel('Atmospheric CO$\mathregular{_{2}}$ concentration (ppm)',
                   fontsize=20, **hfont)
    ax2.xaxis.tick_top()
    ax2.xaxis.set_label_position('top')
    ax2.set_ylim([0, 1.5])
    new_tick_locations = ppm2gtc(np.arange(400, 700, 25))

    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)

    ax2.set_xticklabels(gtc2ppm(new_tick_locations),
                        fontsize=20, **hfont)
    '''
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(ax1.get_xticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    # ax2.xaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(1000, linewidth=2.5, color='k')
    ax1.axhline(0, linewidth=2.5, color='k')
    ax1.axhline(1.5, linewidth=2.5, color='k')

    pink_patch = mpatches.Patch(color='#F5D9CE', label='RCP range')
    red_line = mlines.Line2D([], [], color='#FF0400', linestyle='--',
                             linewidth=2, label='Carbon budget bound')
    black_line = mlines.Line2D([], [], color='k', linestyle='--', linewidth=2,
                               dashes=(10, 10), label='1.5$^{\circ}$C warming')
    b_spr = mpatches.Patch(color='#0708EC', alpha=0.25, linewidth=0)
    ax1.legend(handles=[pink_patch, black_line, mns, (cns, b_spr), red_line,
                        (ubp, spr)],
               labels=['RCP range', '1.5$^{\circ}$C warming', 'RCP2.6 2090s',
                       'Constant WBGT95p',
                       'Carbon budget bound', 'New upper bound'],
               handlelength=3,
               loc=3, scatterpoints=1, prop={'family': 'Arial'},
               bbox_to_anchor=(.75, .78, 2, 1), frameon=True)
    plt.title(title, y=1.01, fontsize=20, **hfont)
    # plt.subplots_adjust(hspace=0, wspace=0.05, top=.95, bottom=0.1, left=.05,
    #                     right=.97)
    un = oub - lb
    p_ch = np.array([float(ub)-oub, float(ub)-float(ubl),
                     float(ubu)-float(ub)])*100/un

    return oub, float(ub), float(ubl), float(ubu), p_ch


def fig4_giorgi(tmeans, conversion, co2plume, rcpplume, coefs):
    from matplotlib import pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.lines as mlines

    def gtc2ppm(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        return ["%.0f" % z for z in x2]

    def gtc2ppmfloat(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        return x2

    def ppm2gtc(x1):
        x2 = np.interp(x1, conversion[:, 1], conversion[:, 2])
        return x2

    def df2gtc(x1):
        dco2 = 390.4*np.exp(np.log(2)*x1/3.74)
        x2 = np.interp(dco2, conversion[:, 1], conversion[:, 2])
        return x2

    def gtc2df(x1):
        x2 = np.interp(x1, conversion[:, 2], conversion[:, 1])
        x3 = (3.74/np.log(2))*np.log(x2/390.4)
        return x3

    def dfdc(c):
        ppm1 = gtc2ppmfloat(c-10)
        ppm2 = gtc2ppmfloat(c+10)
        f1 = 3.74 * np.log(ppm1/278) / np.log(2)
        f2 = 3.74 * np.log(ppm2/278) / np.log(2)
        dfdc = (f2-f1)/20
        return dfdc
    co2 = np.array([390.4, 423.4, 395.8, 550.0])
    tmean = np.zeros((4))
    tmean[0] = np.mean(tmeans['batch_518'])
    tmean[1] = np.mean(tmeans['batch_520'])
    tmean[2] = np.mean(tmeans['batch_521'])
    tmean[3] = np.mean(tmeans['batch_522'])

    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv['d_f'][1] +
                       (x*alpha/beta**2)**2*cv['d_t'][0] -
                       2*(x**2*alpha*cv['d_t'][1]**2/beta**3))
        return sigT

    co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301
    hfont = {'fontname': 'Arial'}
    fig, axs = plt.subplots(3, 7, sharex='col', sharey='row', facecolor='w', edgecolor='k', linewidth=2)
    for i, item in enumerate(coefs):
        alpha = coefs[item][0]
        beta = coefs[item][1]
        cvs = coefs[item][2]
        ax1 = axs[np.floor(i/7), np.remainder(i+7, 7)]
        ax1.set_ylim([0, 1.5])
        ax1.set_xticks(np.arange(0, 1200, 200))
        ax1.set_xlim([0, 650])

        x11 = rcpplume[0, 20]
        y11 = rcpplume[1, 20]
        x22 = rcpplume[0, 3]
        y22 = rcpplume[1, 3]
        g1 = (rcpplume[1, 19]-rcpplume[1, 20])/(rcpplume[0, 19]-rcpplume[0, 20])
        g2 = (rcpplume[1, 4]-rcpplume[1, 3])/(rcpplume[0, 4]-rcpplume[0, 3])
        ax1.axhline(tmean[1]-tmean[0], linestyle='--', color='k', linewidth=2,
                    dashes=(10, 10))
        # ax1.fill_between(co2plume[0, :], co2plume[1, :], facecolor='gray',
        #                  alpha=0.2)
        ax1.fill_between(rcpplume[0, :], rcpplume[1, :], facecolor='#F5D9CE',
                         edgecolor='#eebdaa', alpha=1)
        x = np.arange(0, 1010, .01)
        y = (-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0])
        cns, = ax1.plot(x, y, linewidth=2, color='#0708EC')
        y1 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) -
              1.6445*sigt(gtc2df(x), cvs, alpha, beta))
        y2 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) +
              1.6445*sigt(gtc2df(x), cvs, alpha, beta))
        ax1.fill_between(x, y1, y2, facecolor='#0708EC', edgecolor='#0708EC',
                         alpha=.25)

        oub = (tmean[1]-tmean[0] - y22) / g2 + x22
        lb = (tmean[1]-tmean[0] - y11) / g1 + x11
        ub = x[np.argwhere(np.diff(np.sign(y -
                                   g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu = x[np.argwhere(np.diff(np.sign(y2 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl = x[np.argwhere(np.diff(np.sign(y1 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        spr = ax1.axvspan(ubl, ubu, color='#ff6600', alpha=0.25)
        ax1.axvline(oub, linestyle='--', color='#FF0400', linewidth=2)
        ax1.axvline(lb, linestyle='--', color='#FF0400', linewidth=2)
        ubp = ax1.axvline(ub, linestyle='--', color='#FF7F2A', linewidth=2)
        mns = ax1.scatter(rcp_conc, tmean[1]-tmean[0], color='#5EB3FA',
                          marker='o', s=100, lw=3, label='RCP2.6 2090s')
        '''
        ax2 = ax1.twiny()
        ax2.spines['left'].set_position(('axes', 0.0))
        ax2.spines['right'].set_color('none')
        ax2.spines['bottom'].set_position(('axes', 1.0))
        ax2.spines['top'].set_color('none')
        ax2.spines['bottom'].set_smart_bounds(True)
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.yaxis.set_label_position('left')

        ax2.set_xlabel('Atmospheric CO$\mathregular{_{2}}$ concentration (ppm)',
                       fontsize=20, **hfont)
        ax2.xaxis.tick_top()
        ax2.xaxis.set_label_position('top')
        ax2.set_ylim([0, 1.5])
        new_tick_locations = ppm2gtc(np.arange(400, 700, 25))

        ax2.set_xlim(ax1.get_xlim())
        ax2.set_xticks(new_tick_locations)

        ax2.set_xticklabels(gtc2ppm(new_tick_locations),
                            fontsize=20, **hfont)
        '''
        ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
        ax1.set_xticklabels(ax1.get_xticks(), fontsize=20, **hfont)
        ax1.xaxis.set_tick_params(width=2, length=7.5)
        ax1.yaxis.set_tick_params(width=2, length=7.5)
        # ax2.xaxis.set_tick_params(width=2, length=7.5)
        ax1.axvline(0, linewidth=2.5, color='k')
        ax1.axvline(1000, linewidth=2.5, color='k')
        ax1.axhline(0, linewidth=2.5, color='k')
        ax1.axhline(1.5, linewidth=2.5, color='k')
        '''
        pink_patch = mpatches.Patch(color='#F5D9CE', label='RCP range')
        red_line = mlines.Line2D([], [], color='#FF0400', linestyle='--',
                                 linewidth=2, label='Carbon budget bound')
        black_line = mlines.Line2D([], [], color='k', linestyle='--', linewidth=2,
                                   dashes=(10, 10), label='1.5$^{\circ}$C warming')
        b_spr = mpatches.Patch(color='#0708EC', alpha=0.25, linewidth=0)
        ax1.legend(handles=[pink_patch, black_line, mns, (cns, b_spr), red_line,
                            (ubp, spr)],
                   labels=['RCP range', '1.5$^{\circ}$C warming', 'RCP2.6 2090s',
                           'Constant WBGT95p',
                           'Carbon budget bound', 'New upper bound'],
                   handlelength=3,
                   loc=3, scatterpoints=1, prop={'family': 'Arial'},
                   bbox_to_anchor=(.75, .78, 2, 1), frameon=True)
        '''
        ax1.set_title(item, y=1.01, fontsize=20, **hfont)
    plt.subplots_adjust(hspace=.2, wspace=0.2, top=.95, bottom=0.05, left=.05,
                        right=.95)
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
    plt.ylabel('Global mean temperature anomaly relative to 2006-2016 \
($^{\circ}$C)', fontsize=20, **hfont)
    plt.xlabel('Cumulative total anthropogenic CO$\mathregular{_{2}}$ \
emissions from 2011 (Gt C)', fontsize=20, **hfont)


def boxplot(globprcnts, prcnts):
    from matplotlib import pyplot as plt
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([0.5, 1, 1.5, 3, 3.5, 4, 5.5, 6, 6.5, 9, 9.5, 10,
                  11.5, 12, 12.5, 15, 15.5, 16, 17.5, 18, 18.5])
    globprcnts2 = np.copy(globprcnts)
    globprcnts2[1:, 0] -= globprcnts[0, 0]
    globprcnts2 = np.vstack((globprcnts2[2], globprcnts2[1], globprcnts2[3]))
    prcnts2 = np.copy(prcnts)
    prcnts2[1:, :, 0] -= prcnts[0, :, 0]
    prcnts2 = prcnts2[1:, :]
    hold = np.copy(prcnts2[0, :])
    prcnts2[0, :] = prcnts2[1, :]
    prcnts2[1, :] = hold
    labels = ['', 'Tglob', '', '', 'Tmean', '', '', 'Tmean_trop', '', '',
              'TX90p', '', '', 'TX90p_trop', '', '', 'WBGT95p', '', '',
              'WBGT95p_trop', '']

    y = np.concatenate((globprcnts2[:, 0].flatten(order='F'),
                        prcnts2[:, :, 0].flatten(order='F')))
    yerrs = np.array([np.concatenate((globprcnts2[:, 1].flatten(order='F'),
                                      prcnts2[:, :, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2[:, 2].flatten(order='F'),
                                      prcnts2[:, :, 2].flatten(order='F')))])
    color = ['b', 'k', 'r']
    pp = []
    for i in range(3):
        p = ax1.errorbar(x[i:9:3], y[i:9:3], yerr=yerrs[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(p[2])
    ax2 = ax1.twinx()
    for i in range(3):
        ax2.errorbar(x[9+i:15:3], y[9+i:15:3], yerr=yerrs[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    for i in range(3):
        ax3.errorbar(x[15+i::3], y[15+i::3], yerr=yerrs[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)

    ax2.spines['left'].set_position(('data', 8.5))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 14.5))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax1.set_xlim([0, 19])
    ax1.set_ylabel('Temperature anomaly [K]', fontweight='bold', fontsize=20)
    ax2.set_ylabel('Percentage anomaly [%]', fontweight='bold', fontsize=20)
    ax3.set_ylabel('Degree month anomaly [K months]', fontweight='bold',
                   fontsize=20)

    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontweight='bold', fontsize=20)

    ax1.set_xticklabels(labels, fontweight='bold', fontsize=20,
                        rotation='vertical')
    ax2.set_yticklabels(ax2.get_yticks(), fontweight='bold', fontsize=20)
    ax3.set_yticklabels(ax3.get_yticks(), fontweight='bold', fontsize=20)

    ax1.legend(handles=pp, labels=['Low CO2', 'Mean CO2', 'High CO2'],
               loc='upper right', scatterpoints=1)

    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.2, left=.05,
                        right=.99)


def boxplotprecip(globprcnts, prcnts):
    from matplotlib import pyplot as plt
    hfont = {'fontname': 'Arial'}
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([0.5, 1, 1.5, 3, 3.5, 4, 5.5, 6, 6.5, 9, 9.5, 10,
                  11.5, 12, 12.5, 15, 15.5, 16, 17.5, 18, 18.5, 21, 21.5, 22,
                  23.5, 24, 24.5, 26, 26.5, 27, 29.5, 30, 30.5, 32, 32.5, 33])
    globprcnts2 = np.copy(globprcnts)
    globprcnts2[0, 1:, 0] -= globprcnts[0, 0, 0]
    globprcnts2 = np.vstack((globprcnts2[0, 2], globprcnts2[0, 1],
                             globprcnts2[0, 3]))
    globprcnts3 = np.copy(globprcnts)
    globprcnts3[1, 1:, 0] -= globprcnts[1, 0, 0]
    globprcnts3 = np.vstack((globprcnts3[1, 2], globprcnts3[1, 1],
                             globprcnts3[1, 3]))
    prcnts2 = np.copy(prcnts)
    prcnts2[1:, :, 0] -= prcnts[0, :, 0]
    prcnts2 = prcnts2[1:, :]
    hold = np.copy(prcnts2[0, :])
    prcnts2[0, :] = prcnts2[1, :]
    prcnts2[1, :] = hold
    labels = ['', '$\overline{T}_{global}$', '', '', 'Tmean_NH', '', '', 'Tmean_trop', '', '',
              'TX90p_NH', '', '', 'TX90p_trop', '', '', 'WBGT95p_NH', '', '',
              'WBGT95p_trop', '', '', 'Rglob', '', '', 'Rmean_NH', '', '',
              'Rmean_trop',  '', '', 'R95p_NH', '', '', 'R95p_trop', '']

    y = np.concatenate((globprcnts2[:, 0].flatten(order='F'),
                        prcnts2[:, :6, 0].flatten(order='F')))
    yerrs = np.array([np.concatenate((globprcnts2[:, 1].flatten(order='F'),
                                      prcnts2[:, :6, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2[:, 2].flatten(order='F'),
                                      prcnts2[:, :6, 2].flatten(order='F')))])
    y1 = np.concatenate((globprcnts3[:, 0].flatten(order='F'),
                        prcnts2[:, 6:, 0].flatten(order='F')))
    yerrs1 = np.array([np.concatenate((globprcnts3[:, 1].flatten(order='F'),
                                      prcnts2[:, 6:, 1].flatten(order='F'))),
                      np.concatenate((globprcnts3[:, 2].flatten(order='F'),
                                      prcnts2[:, 6:, 2].flatten(order='F')))])
    color = ['b', 'k', 'r']
    pp = []
    for i in range(3):
        p = ax1.errorbar(x[i:9:3], y[i:9:3], yerr=yerrs[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(p[2])
    ax2 = ax1.twinx()
    for i in range(3):
        ax2.errorbar(x[9+i:15:3], y[9+i:15:3], yerr=yerrs[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    for i in range(3):
        ax3.errorbar(x[15+i:21:3], y[15+i::3], yerr=yerrs[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax4 = ax1.twinx()
    for i in range(3):
        ax4.errorbar(x[21+i:30:3], y1[i:9:3]*360, yerr=yerrs1[:, i:9:3]*360,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax5 = ax1.twinx()
    for i in range(3):
        ax5.errorbar(x[30+i::3], y1[9+i::3], yerr=yerrs1[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)

    ax2.spines['left'].set_position(('data', 8.5))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 14.5))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax4.spines['left'].set_position(('data', 20.5))
    ax4.spines['right'].set_color('none')
    ax4.spines['bottom'].set_position(('axes', 0.0))
    ax4.spines['top'].set_color('none')
    ax4.spines['bottom'].set_smart_bounds(True)
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')

    ax5.spines['left'].set_position(('data', 29))
    ax5.spines['right'].set_color('none')
    ax5.spines['bottom'].set_position(('axes', 0.0))
    ax5.spines['top'].set_color('none')
    ax5.spines['bottom'].set_smart_bounds(True)
    ax5.xaxis.set_ticks_position('bottom')
    ax5.yaxis.set_ticks_position('left')
    ax5.yaxis.set_label_position('left')

    ax1.set_xlim([0, 33.5])
    ax1.set_ylim([.6, 1.8])
    ax2.set_ylim([3.75, 14.75])
    ax3.set_ylim([.175, .675])
    ax4.set_ylim([-15, 45])
    ax5.set_ylim([-0.025, 0.525])
    ax1.set_ylabel('Temperature anomaly ($^{\circ}$C)', fontsize=20, **hfont)
    ax2.set_ylabel('Anomaly (%)', fontsize=20,
                   **hfont)
    ax3.set_ylabel('Degree month anomaly ($^{\circ}$C months)',
                   fontsize=20, **hfont)
    ax4.set_ylabel('Precipitation anomaly (mm yr$\mathregular{^{-1}}$)',
                   fontsize=20, **hfont)
    ax5.set_ylabel('Anomaly (%)',
                   fontsize=20, **hfont)

    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation='vertical', **hfont)
    ax2.set_yticklabels(ax2.get_yticks(), fontsize=20, **hfont)
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=20, **hfont)
    ax4.set_yticklabels(ax4.get_yticks(), fontsize=20, **hfont)
    ax5.set_yticklabels(ax5.get_yticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax2.yaxis.set_tick_params(width=2, length=7.5)
    ax3.yaxis.set_tick_params(width=2, length=7.5)
    ax4.yaxis.set_tick_params(width=2, length=7.5)
    ax5.yaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(8.5, linewidth=1.5, color='k')
    ax1.axvline(14.5, linewidth=1.5, color='k')
    ax1.axvline(20.5, linewidth=1.5, color='k')
    ax1.axvline(29, linewidth=1.5, color='k')
    ax1.axvline(33.5, linewidth=2.5, color='k')
    ax1.axhline(0.6, linewidth=2.5, color='k')
    ax1.axhline(1.8, linewidth=2.5, color='k')
    ax1.legend(handles=pp, labels=['Low CO$\mathregular{_{2}}$',
                                   'Mean CO$\mathregular{_{2}}$',
                                   'High CO$\mathregular{_{2}}$'],
               loc=3, scatterpoints=1, prop={'family': 'Arial'},
               bbox_to_anchor=(.885, .86, 1., .102))
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.2, left=.05,
                        right=.99)


def boxplotprecip_days(globprcnts, prcnts):
    from matplotlib import pyplot as plt
    hfont = {'fontname': 'Arial'}
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([0.5, 1, 1.5, 3, 3.5, 4, 5.5, 6, 6.5, 9, 9.5, 10,
                  11.5, 12, 12.5, 15, 15.5, 16, 17.5, 18, 18.5, 21, 21.5, 22,
                  23.5, 24, 24.5, 26, 26.5, 27, 29.5, 30, 30.5, 32, 32.5, 33])
    globprcnts2 = np.copy(globprcnts)
    globprcnts2[0, 1:, 0] -= globprcnts[0, 0, 0]
    globprcnts2 = np.vstack((globprcnts2[0, 2], globprcnts2[0, 1],
                             globprcnts2[0, 3]))
    globprcnts3 = np.copy(globprcnts)
    globprcnts3[2, 1:, 0] -= globprcnts[2, 0, 0]
    globprcnts3 = np.vstack((globprcnts3[2, 2], globprcnts3[2, 1],
                             globprcnts3[2, 3]))
    prcnts2 = np.copy(prcnts)
    prcnts2[1:, :, 0] -= prcnts[0, :, 0]
    prcnts2 = prcnts2[1:, :]
    hold = np.copy(prcnts2[0, :])
    prcnts2[0, :] = prcnts2[1, :]
    prcnts2[1, :] = hold
    labels = ['', '$\mathregular{\overline{T}_{GLOB}}$', '', '', '$\mathregular{\overline{T}_{NH}}$', '', '', '$\mathregular{\overline{T}_{TROP}}$', '', '',
              '$\mathregular{\overline{TX90p}_{NH}}$', '', '', '$\mathregular{\overline{TX90p}_{TROP}}$', '', '', '$\mathregular{\overline{WBGT95p}_{NH}}$', '', '',
              '$\mathregular{\overline{WBGT95p}_{TROP}}$', '', '', '$\mathregular{\overline{R}_{GLOB}}$', '', '', '$\mathregular{\overline{R}_{NH}}$', '', '',
              '$\mathregular{\overline{R}_{TROP}}$',  '', '', '$\mathregular{\overline{R95p}_{NH}}$', '', '', '$\mathregular{\overline{R95p}_{TROP}}$', '']

    y = np.concatenate((globprcnts2[:, 0].flatten(order='F'),
                        prcnts2[:, :6, 0].flatten(order='F')))
    yerrs = np.array([np.concatenate((globprcnts2[:, 1].flatten(order='F'),
                                      prcnts2[:, :6, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2[:, 2].flatten(order='F'),
                                      prcnts2[:, :6, 2].flatten(order='F')))])
    y1 = np.concatenate((globprcnts3[:, 0].flatten(order='F'),
                        prcnts2[:, 6:, 0].flatten(order='F')))
    yerrs1 = np.array([np.concatenate((globprcnts3[:, 1].flatten(order='F'),
                                      prcnts2[:, 6:, 1].flatten(order='F'))),
                      np.concatenate((globprcnts3[:, 2].flatten(order='F'),
                                      prcnts2[:, 6:, 2].flatten(order='F')))])
    color = ['b', 'k', 'r']
    pp = []
    for i in range(3):
        p = ax1.errorbar(x[i:9:3], y[i:9:3], yerr=yerrs[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(p[2])
    ax2 = ax1.twinx()
    for i in range(3):
        ax2.errorbar(x[9+i:15:3], y[9+i:15:3], yerr=yerrs[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    for i in range(3):
        ax3.errorbar(x[15+i:21:3], y[15+i::3], yerr=yerrs[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax4 = ax1.twinx()
    for i in range(3):
        ax4.errorbar(x[21+i:30:3], y1[i:9:3]*360, yerr=yerrs1[:, i:9:3]*360,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax5 = ax1.twinx()
    for i in range(3):
        ax5.errorbar(x[30+i::3], y1[9+i::3], yerr=yerrs1[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)

    ax2.spines['left'].set_position(('data', 8.5))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 14.5))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax4.spines['left'].set_position(('data', 20.5))
    ax4.spines['right'].set_color('none')
    ax4.spines['bottom'].set_position(('axes', 0.0))
    ax4.spines['top'].set_color('none')
    ax4.spines['bottom'].set_smart_bounds(True)
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')

    ax5.spines['left'].set_position(('data', 29))
    ax5.spines['right'].set_color('none')
    ax5.spines['bottom'].set_position(('axes', 0.0))
    ax5.spines['top'].set_color('none')
    ax5.spines['bottom'].set_smart_bounds(True)
    ax5.xaxis.set_ticks_position('bottom')
    ax5.yaxis.set_ticks_position('left')
    ax5.yaxis.set_label_position('left')

    ax1.set_xlim([0, 33.5])
    ax1.set_ylim([.6, 1.8])
    ax2.set_ylim([3.75, 13.25])
    ax3.set_ylim([.175, .675])
    ax4.set_ylim([-15, 45])
    ax5.set_ylim([-0.025, .475])
    ax1.set_ylabel('Temperature anomaly ($^{\circ}$C)', fontsize=20, **hfont)
    ax2.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)', fontsize=20,
                   **hfont)
    ax3.set_ylabel('Anomaly ($^{\circ}$C months)',
                   fontsize=20, **hfont)
    ax4.set_ylabel('Precipitation anomaly (mm year$\mathregular{^{-1}}$)',
                   fontsize=20, **hfont)
    ax5.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)',
                   fontsize=20, **hfont)

    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation='vertical', **hfont)
    ax2.set_yticklabels(ax2.get_yticks(), fontsize=20, **hfont)
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=20, **hfont)
    ax4.set_yticklabels(ax4.get_yticks(), fontsize=20, **hfont)
    ax5.set_yticklabels(ax5.get_yticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax2.yaxis.set_tick_params(width=2, length=7.5)
    ax3.yaxis.set_tick_params(width=2, length=7.5)
    ax4.yaxis.set_tick_params(width=2, length=7.5)
    ax5.yaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(8.5, linewidth=1.5, color='k')
    ax1.axvline(14.5, linewidth=1.5, color='k')
    ax1.axvline(20.5, linewidth=1.5, color='k')
    ax1.axvline(29, linewidth=1.5, color='k')
    ax1.axvline(33.5, linewidth=2.5, color='k')
    ax1.axhline(0.6, linewidth=2.5, color='k')
    ax1.axhline(1.8, linewidth=2.5, color='k')
    ax1.legend(handles=pp, labels=['Low CO$\mathregular{_{2}}$',
                                   'Mean CO$\mathregular{_{2}}$',
                                   'High CO$\mathregular{_{2}}$'],
               loc=3, scatterpoints=1, prop={'family': 'Arial'},
               bbox_to_anchor=(.885, .86, 1., .102))
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.2, left=.05,
                        right=.99)


def plotall(data, colorlimit, meaning='mean', mask='yes', precip='no',
            cbarleg='cbarlg', pltlbl='pltlbl', b1='batch_521', b2='batch_522'):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')
    if meaning == 'mean':
        winter = np.mean(data[b2][0, :], axis=0) - np.mean(data[b1][0, :], axis=0)
        summer = np.mean(data[b2][1, :], axis=0) - np.mean(data[b1][1, :], axis=0)
        sig0 = ttest(data[b1][0, :],
                     data[b2][0, :], siglevel=10)
        sig1 = ttest(data[b1][1, :],
                     data[b2][1, :], siglevel=10)
    elif meaning == 'non':
        winter = data[b2][0, :] - data[b1][0, :]
        summer = data[b2][1, :] - data[b1][1, :]
        sig0 = np.zeros((145, 192))
        sig1 = np.zeros((145, 192))
    else:
        winter = np.mean(data[b2][meaning][0, :], axis=0) - np.mean(data[b1][meaning][0, :], axis=0)
        summer = np.mean(data[b2][meaning][1, :], axis=0) - np.mean(data[b1][meaning][1, :], axis=0)
        sig0 = ttest(data[b2][meaning][0, :],
                     data[b1][meaning][0, :], siglevel=10)
        sig1 = ttest(data[b2][meaning][1, :],
                     data[b1][meaning][1, :], siglevel=10)
    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        winter = np.ma.masked_array(winter, mask=np.logical_not(lsm))
        summer = np.ma.masked_array(summer, mask=np.logical_not(lsm))
        sig0 = np.ma.masked_array(sig0, mask=np.logical_not(lsm))
        sig1 = np.ma.masked_array(sig1, mask=np.logical_not(lsm))

    winter, lon1 = shiftgrid(180., winter, lon, start=False)
    winter, lon1 = addcyclic(winter, lon1)
    sig0, lon2 = shiftgrid(180., sig0, lon, start=False)
    sig0, lon2 = addcyclic(sig0, lon2)
    summer, lon3 = shiftgrid(180., summer, lon, start=False)
    summer, lon3 = addcyclic(summer, lon3)
    sig1, lon4 = shiftgrid(180., sig1, lon, start=False)
    sig1, lon4 = addcyclic(sig1, lon4)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrs = np.linspace(-colorlimit, colorlimit, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    if precip == 'yes':
        my_cmap = my_cmap[::-1]
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap

    ax1 = fig.add_subplot(1, 2, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, winter, ctrs,
                      cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    plt.contourf(x, y, sig0, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax2 = fig.add_subplot(1, 2, 2)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot1 = m.contourf(x, y, summer, ctrs,
                       cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                       extend='both')
    # parallels = m.drawparallels(np.arange(-90., 91., 15.))
    # meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    # m.drawparallels(parallels, labels=[True, True, True, True])
    # m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.contourf(x, y, sig1, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    hfont = {'fontname': 'Arial'}
    cbar_ax = fig.add_axes([0.2, 0.3, 0.6, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max')
    b.set_label(label=cbarleg, size=20, fontsize=20, fontname='Arial')
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)

    pad = 15  # in points
    ax1.annotate('DJF', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('JJA', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax1.annotate(pltlbl, xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.95)


def fig1(datat, dataw, datar):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    winterr = np.mean(datar['batch_522'][0, :], axis=0) - np.mean(datar['batch_521'][0, :], axis=0)
    summerr = np.mean(datar['batch_522'][1, :], axis=0) - np.mean(datar['batch_521'][1, :], axis=0)
    sig0r = ttest(datar['batch_521'][0, :],
                  datar['batch_522'][0, :], siglevel=10)
    sig1r = ttest(datar['batch_521'][1, :],
                  datar['batch_522'][1, :], siglevel=10)

    wintert = np.mean(datat['batch_522'][2][0, :], axis=0) - np.mean(datat['batch_521'][2][0, :], axis=0)
    summert = np.mean(datat['batch_522'][2][1, :], axis=0) - np.mean(datat['batch_521'][2][1, :], axis=0)
    sig0t = ttest(datat['batch_522'][2][0, :],
                  datat['batch_521'][2][0, :], siglevel=10)
    sig1t = ttest(datat['batch_522'][2][1, :],
                  datat['batch_521'][2][1, :], siglevel=10)

    winterw = np.mean(dataw['batch_522'][0, :], axis=0) - np.mean(dataw['batch_521'][0, :], axis=0)
    summerw = np.mean(dataw['batch_522'][1, :], axis=0) - np.mean(dataw['batch_521'][1, :], axis=0)
    sig0w = ttest(dataw['batch_522'][0, :],
                  dataw['batch_521'][0, :], siglevel=10)
    sig1w = ttest(dataw['batch_522'][1, :],
                  dataw['batch_521'][1, :], siglevel=10)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    wintert = np.ma.masked_array(wintert, mask=np.logical_not(lsm))
    summert = np.ma.masked_array(summert, mask=np.logical_not(lsm))
    sig0t = np.ma.masked_array(sig0t, mask=np.logical_not(lsm))
    sig1t = np.ma.masked_array(sig1t, mask=np.logical_not(lsm))
    winterw = np.ma.masked_array(winterw, mask=np.logical_not(lsm))
    summerw = np.ma.masked_array(summerw, mask=np.logical_not(lsm))
    sig0w = np.ma.masked_array(sig0w, mask=np.logical_not(lsm))
    sig1w = np.ma.masked_array(sig1w, mask=np.logical_not(lsm))

    wintert, lon1 = shiftgrid(180., wintert, lon, start=False)
    wintert, lon1 = addcyclic(wintert, lon1)
    sig0t, lon2 = shiftgrid(180., sig0t, lon, start=False)
    sig0t, lon2 = addcyclic(sig0t, lon2)
    summert, lon3 = shiftgrid(180., summert, lon, start=False)
    summert, lon3 = addcyclic(summert, lon3)
    sig1t, lon4 = shiftgrid(180., sig1t, lon, start=False)
    sig1t, lon4 = addcyclic(sig1t, lon4)
    winterr, lon1 = shiftgrid(180., winterr, lon, start=False)
    winterr, lon1 = addcyclic(winterr, lon1)
    sig0r, lon2 = shiftgrid(180., sig0r, lon, start=False)
    sig0r, lon2 = addcyclic(sig0r, lon2)
    summerr, lon3 = shiftgrid(180., summerr, lon, start=False)
    summerr, lon3 = addcyclic(summerr, lon3)
    sig1r, lon4 = shiftgrid(180., sig1r, lon, start=False)
    sig1r, lon4 = addcyclic(sig1r, lon4)
    winterw, lon1 = shiftgrid(180., winterw, lon, start=False)
    winterw, lon1 = addcyclic(winterw, lon1)
    sig0w, lon2 = shiftgrid(180., sig0w, lon, start=False)
    sig0w, lon2 = addcyclic(sig0w, lon2)
    summerw, lon3 = shiftgrid(180., summerw, lon, start=False)
    summerw, lon3 = addcyclic(summerw, lon3)
    sig1w, lon4 = shiftgrid(180., sig1w, lon, start=False)
    sig1w, lon4 = addcyclic(sig1w, lon4)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrst = np.linspace(-15, 15, 17)
    ctrsw = np.linspace(-1, 1, 17)
    ctrsr = np.linspace(-2, 2, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap1 = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap1)
    cmp1 = newcmap1

    ax1 = fig.add_subplot(3, 2, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, wintert, ctrst,
                      cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                      extend='both')
    plt.contourf(x, y, sig0t, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20, y=0.9)
    ax2 = fig.add_subplot(3, 2, 2)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot2 = m.contourf(x, y, summert, ctrst,
                       cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                       extend='both')
    plt.contourf(x, y, sig1t, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20, y=0.9)
    ax3 = fig.add_subplot(3, 2, 3)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot3 = m.contourf(x, y, winterw, ctrsw,
                       cmap=cmp, vmin=np.min(ctrsw), vmax=np.max(ctrsw),
                       extend='both')
    plt.contourf(x, y, sig0w, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20, y=0.9)
    ax4 = fig.add_subplot(3, 2, 4)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot4 = m.contourf(x, y, summerw, ctrsw,
                       cmap=cmp, vmin=np.min(ctrsw), vmax=np.max(ctrsw),
                       extend='both')
    plt.contourf(x, y, sig1w, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20, y=0.9)
    ax5 = fig.add_subplot(3, 2, 5)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot5 = m.contourf(x, y, winterr, ctrsr,
                       cmap=cmp1, vmin=np.min(ctrsr), vmax=np.max(ctrsr),
                       extend='both')
    plt.contourf(x, y, sig0r, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax5.set_title('e', loc='left', fontname='Arial', fontsize=20, y=0.9)
    ax6 = fig.add_subplot(3, 2, 6)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot6 = m.contourf(x, y, summerr, ctrsr,
                       cmap=cmp1, vmin=np.min(ctrsr), vmax=np.max(ctrsr),
                       extend='both')
    plt.contourf(x, y, sig1r, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax6.set_title('f', loc='left', fontname='Arial', fontsize=20, y=0.9)

    cbar_ax = fig.add_axes([0.925, 0.7125, 0.005, 0.2])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max')
    b.set_label(label='Difference (%)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.405, 0.005, 0.2])
    b = fig.colorbar(plot4, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max')
    b.set_label(label='Difference ($^{\circ}$C months)', size=20,
                fontsize=20, fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.1, 0.005, 0.2])
    b = fig.colorbar(plot6, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max')
    b.set_label(label='Difference (%)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 15  # in points
    ax1.annotate('DJF', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('JJA', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax1.annotate('TX90p', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('WBGT95p', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax5.annotate('R95p', xy=(0, 0.5), xytext=(-ax5.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax5.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0.05, wspace=0, top=.95, bottom=0.05,
                        left=.05, right=.9)


def fig1_days(datat, dataw, datar):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    winterr = np.mean(datar['batch_522'][0, :], axis=0) - np.mean(datar['batch_521'][0, :], axis=0)
    summerr = np.mean(datar['batch_522'][1, :], axis=0) - np.mean(datar['batch_521'][1, :], axis=0)
    sig0r = ttest(datar['batch_521'][0, :],
                  datar['batch_522'][0, :], siglevel=10)
    sig1r = ttest(datar['batch_521'][1, :],
                  datar['batch_522'][1, :], siglevel=10)

    wintert = np.mean(datat['batch_522'][0, :], axis=0) - np.mean(datat['batch_521'][0, :], axis=0)
    summert = np.mean(datat['batch_522'][1, :], axis=0) - np.mean(datat['batch_521'][1, :], axis=0)
    sig0t = ttest(datat['batch_522'][0, :],
                  datat['batch_521'][0, :], siglevel=10)
    sig1t = ttest(datat['batch_522'][1, :],
                  datat['batch_521'][1, :], siglevel=10)

    winterw = np.mean(dataw['batch_522'][0, :], axis=0) - np.mean(dataw['batch_521'][0, :], axis=0)
    summerw = np.mean(dataw['batch_522'][1, :], axis=0) - np.mean(dataw['batch_521'][1, :], axis=0)
    sig0w = ttest(dataw['batch_522'][0, :],
                  dataw['batch_521'][0, :], siglevel=10)
    sig1w = ttest(dataw['batch_522'][1, :],
                  dataw['batch_521'][1, :], siglevel=10)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    wintert = np.ma.masked_array(wintert, mask=np.logical_not(lsm))
    summert = np.ma.masked_array(summert, mask=np.logical_not(lsm))
    sig0t = np.ma.masked_array(sig0t, mask=np.logical_not(lsm))
    sig1t = np.ma.masked_array(sig1t, mask=np.logical_not(lsm))
    winterw = np.ma.masked_array(winterw, mask=np.logical_not(lsm))
    summerw = np.ma.masked_array(summerw, mask=np.logical_not(lsm))
    sig0w = np.ma.masked_array(sig0w, mask=np.logical_not(lsm))
    sig1w = np.ma.masked_array(sig1w, mask=np.logical_not(lsm))

    wintert, lon1 = shiftgrid(180., wintert, lon, start=False)
    wintert, lon1 = addcyclic(wintert, lon1)
    sig0t, lon2 = shiftgrid(180., sig0t, lon, start=False)
    sig0t, lon2 = addcyclic(sig0t, lon2)
    summert, lon3 = shiftgrid(180., summert, lon, start=False)
    summert, lon3 = addcyclic(summert, lon3)
    sig1t, lon4 = shiftgrid(180., sig1t, lon, start=False)
    sig1t, lon4 = addcyclic(sig1t, lon4)
    winterr, lon1 = shiftgrid(180., winterr, lon, start=False)
    winterr, lon1 = addcyclic(winterr, lon1)
    sig0r, lon2 = shiftgrid(180., sig0r, lon, start=False)
    sig0r, lon2 = addcyclic(sig0r, lon2)
    summerr, lon3 = shiftgrid(180., summerr, lon, start=False)
    summerr, lon3 = addcyclic(summerr, lon3)
    sig1r, lon4 = shiftgrid(180., sig1r, lon, start=False)
    sig1r, lon4 = addcyclic(sig1r, lon4)
    winterw, lon1 = shiftgrid(180., winterw, lon, start=False)
    winterw, lon1 = addcyclic(winterw, lon1)
    sig0w, lon2 = shiftgrid(180., sig0w, lon, start=False)
    sig0w, lon2 = addcyclic(sig0w, lon2)
    summerw, lon3 = shiftgrid(180., summerw, lon, start=False)
    summerw, lon3 = addcyclic(summerw, lon3)
    sig1w, lon4 = shiftgrid(180., sig1w, lon, start=False)
    sig1w, lon4 = addcyclic(sig1w, lon4)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrst = np.linspace(-15, 15, 17)
    ctrsw = np.linspace(-1, 1, 17)
    ctrsr = np.linspace(-2, 2, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap1 = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap1)
    cmp1 = newcmap1

    ax1 = fig.add_subplot(3, 2, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, wintert, ctrst,
                      cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                      extend='both')
    plt.contourf(x, y, sig0t, levels=[0, .5, 1.5],
                 hatches=["", "."], alpha=0)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax2 = fig.add_subplot(3, 2, 2)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot2 = m.contourf(x, y, summert, ctrst,
                       cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                       extend='both')
    plt.contourf(x, y, sig1t, levels=[0, .5, 1.5],
                 hatches=["", "."], alpha=0)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax3 = fig.add_subplot(3, 2, 3)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot3 = m.contourf(x, y, winterw, ctrsw,
                       cmap=cmp, vmin=np.min(ctrsw), vmax=np.max(ctrsw),
                       extend='both')
    plt.contourf(x, y, sig0w, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax4 = fig.add_subplot(3, 2, 4)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot4 = m.contourf(x, y, summerw, ctrsw,
                       cmap=cmp, vmin=np.min(ctrsw), vmax=np.max(ctrsw),
                       extend='both')
    plt.contourf(x, y, sig1w, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax5 = fig.add_subplot(3, 2, 5)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot5 = m.contourf(x, y, winterr, ctrsr,
                       cmap=cmp1, vmin=np.min(ctrsr), vmax=np.max(ctrsr),
                       extend='both')
    plt.contourf(x, y, sig0r, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax5.set_title('e', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax6 = fig.add_subplot(3, 2, 6)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot6 = m.contourf(x, y, summerr, ctrsr,
                       cmap=cmp1, vmin=np.min(ctrsr), vmax=np.max(ctrsr),
                       extend='both')
    plt.contourf(x, y, sig1r, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax6.set_title('f', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)

    cbar_ax = fig.add_axes([0.925, 0.7125, 0.005, 0.2])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (days season$\mathregular{^{-1}}$)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.405, 0.005, 0.2])
    b = fig.colorbar(plot4, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference ($^{\circ}$C months)', size=20,
                fontsize=20, fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.1, 0.005, 0.2])
    b = fig.colorbar(plot6, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (days season$\mathregular{^{-1}}$)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 15  # in points
    ax1.annotate('DJF', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('JJA', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax1.annotate('TX90p', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('WBGT95p', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax5.annotate('R95p', xy=(0, 0.5), xytext=(-ax5.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax5.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0.05, wspace=0, top=.95, bottom=0.05,
                        left=.05, right=.9)


def ttest(series1, series2, siglevel=10, testtype='two'):
    """
    Student t-test for a time series

    Parameters
    ----------
    series1: array
        control run
    series2: array
        forced run

    Returns
    -------
    sig: bool
        is it significant
    sig: bool
        is it discernible and significant
    """
    import scipy.stats as st
    sigarray = np.full((np.ma.size(series1, axis=1), 192), siglevel)
    if testtype == 'two':
        a = 2
    elif testtype == 'one':
        a = 1
    else:
        print("Error, test type must be 'one' or 'two'")
    d = np.sqrt(np.var(series1, axis=0, ddof=1) + np.var(series2, axis=0,
                                                         ddof=1))
    z1 = (np.mean(series1, axis=0) - np.mean(series2, axis=0)) / d
    p1 = 1 - st.norm.cdf(np.abs(z1))
    sig = np.greater_equal(sigarray, p1*100*a).astype(int)
    return sig


def fig1_extra(datat, datatx, dataw, datar):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    winterr = np.mean(datar['batch_522'][0, :], axis=0) - np.mean(datar['batch_521'][0, :], axis=0)
    summerr = np.mean(datar['batch_522'][1, :], axis=0) - np.mean(datar['batch_521'][1, :], axis=0)
    sig0r = ttest(datar['batch_521'][0, :],
                  datar['batch_522'][0, :], siglevel=10)
    sig1r = ttest(datar['batch_521'][1, :],
                  datar['batch_522'][1, :], siglevel=10)

    wintert = np.mean(datat['batch_522'][0, :], axis=0) - np.mean(datat['batch_521'][0, :], axis=0)
    summert = np.mean(datat['batch_522'][1, :], axis=0) - np.mean(datat['batch_521'][1, :], axis=0)
    sig0t = ttest(datat['batch_522'][0, :],
                  datat['batch_521'][0, :], siglevel=10)
    sig1t = ttest(datat['batch_522'][1, :],
                  datat['batch_521'][1, :], siglevel=10)

    winterw = np.mean(dataw['batch_522'][0, :], axis=0) - np.mean(dataw['batch_521'][0, :], axis=0)
    summerw = np.mean(dataw['batch_522'][1, :], axis=0) - np.mean(dataw['batch_521'][1, :], axis=0)
    sig0w = ttest(dataw['batch_522'][0, :],
                  dataw['batch_521'][0, :], siglevel=10)
    sig1w = ttest(dataw['batch_522'][1, :],
                  dataw['batch_521'][1, :], siglevel=10)

    wintertx = np.mean(datatx['batch_522'][0][0, :], axis=0) - np.mean(datatx['batch_521'][0][0, :], axis=0)
    summertx = np.mean(datatx['batch_522'][0][1, :], axis=0) - np.mean(datatx['batch_521'][0][1, :], axis=0)
    sig0tx = ttest(datatx['batch_522'][0][0, :],
                  datatx['batch_521'][0][0, :], siglevel=10)
    sig1tx = ttest(datatx['batch_522'][0][1, :],
                  datatx['batch_521'][0][1, :], siglevel=10)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    wintert = np.ma.masked_array(wintert, mask=np.logical_not(lsm))
    summert = np.ma.masked_array(summert, mask=np.logical_not(lsm))
    sig0t = np.ma.masked_array(sig0t, mask=np.logical_not(lsm))
    sig1t = np.ma.masked_array(sig1t, mask=np.logical_not(lsm))
    winterw = np.ma.masked_array(winterw, mask=np.logical_not(lsm))
    summerw = np.ma.masked_array(summerw, mask=np.logical_not(lsm))
    sig0w = np.ma.masked_array(sig0w, mask=np.logical_not(lsm))
    sig1w = np.ma.masked_array(sig1w, mask=np.logical_not(lsm))
    wintertx = np.ma.masked_array(wintertx, mask=np.logical_not(lsm))
    summertx = np.ma.masked_array(summertx, mask=np.logical_not(lsm))
    sig0tx = np.ma.masked_array(sig0tx, mask=np.logical_not(lsm))
    sig1tx = np.ma.masked_array(sig1tx, mask=np.logical_not(lsm))

    wintert, lon1 = shiftgrid(180., wintert, lon, start=False)
    wintert, lon1 = addcyclic(wintert, lon1)
    sig0t, lon2 = shiftgrid(180., sig0t, lon, start=False)
    sig0t, lon2 = addcyclic(sig0t, lon2)
    summert, lon3 = shiftgrid(180., summert, lon, start=False)
    summert, lon3 = addcyclic(summert, lon3)
    sig1t, lon4 = shiftgrid(180., sig1t, lon, start=False)
    sig1t, lon4 = addcyclic(sig1t, lon4)
    winterr, lon1 = shiftgrid(180., winterr, lon, start=False)
    winterr, lon1 = addcyclic(winterr, lon1)
    sig0r, lon2 = shiftgrid(180., sig0r, lon, start=False)
    sig0r, lon2 = addcyclic(sig0r, lon2)
    summerr, lon3 = shiftgrid(180., summerr, lon, start=False)
    summerr, lon3 = addcyclic(summerr, lon3)
    sig1r, lon4 = shiftgrid(180., sig1r, lon, start=False)
    sig1r, lon4 = addcyclic(sig1r, lon4)
    winterw, lon1 = shiftgrid(180., winterw, lon, start=False)
    winterw, lon1 = addcyclic(winterw, lon1)
    sig0w, lon2 = shiftgrid(180., sig0w, lon, start=False)
    sig0w, lon2 = addcyclic(sig0w, lon2)
    summerw, lon3 = shiftgrid(180., summerw, lon, start=False)
    summerw, lon3 = addcyclic(summerw, lon3)
    sig1w, lon4 = shiftgrid(180., sig1w, lon, start=False)
    sig1w, lon4 = addcyclic(sig1w, lon4)
    wintertx, lon1 = shiftgrid(180., wintertx, lon, start=False)
    wintertx, lon1 = addcyclic(wintertx, lon1)
    sig0tx, lon2 = shiftgrid(180., sig0tx, lon, start=False)
    sig0tx, lon2 = addcyclic(sig0tx, lon2)
    summertx, lon3 = shiftgrid(180., summertx, lon, start=False)
    summertx, lon3 = addcyclic(summertx, lon3)
    sig1tx, lon4 = shiftgrid(180., sig1tx, lon, start=False)
    sig1tx, lon4 = addcyclic(sig1tx, lon4)
    meshlon, meshlat = np.meshgrid(lon1, lat)
    ctrst = np.linspace(-2, 2, 17)
    ctrsw = np.linspace(-1, 1, 17)
    ctrsr = np.linspace(-1, 1, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap1 = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap1)
    cmp1 = newcmap1

    ax1 = fig.add_subplot(4, 2, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, wintert, ctrst,
                      cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                      extend='both')
    plt.contourf(x, y, sig0t, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax2 = fig.add_subplot(4, 2, 2)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot2 = m.contourf(x, y, summert, ctrst,
                       cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                       extend='both')
    plt.contourf(x, y, sig1t, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax3 = fig.add_subplot(4, 2, 3)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot3 = m.contourf(x, y, wintertx, ctrst,
                       cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                       extend='both')
    plt.contourf(x, y, sig0tx, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax4 = fig.add_subplot(4, 2, 4)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot4 = m.contourf(x, y, summertx, ctrst,
                       cmap=cmp, vmin=np.min(ctrst), vmax=np.max(ctrst),
                       extend='both')
    plt.contourf(x, y, sig1tx, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax5 = fig.add_subplot(4, 2, 5)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot5 = m.contourf(x, y, winterw, ctrsw,
                       cmap=cmp, vmin=np.min(ctrsw), vmax=np.max(ctrsw),
                       extend='both')
    plt.contourf(x, y, sig0w, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax5.set_title('e', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax6 = fig.add_subplot(4, 2, 6)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot6 = m.contourf(x, y, summerw, ctrsw,
                       cmap=cmp, vmin=np.min(ctrsw), vmax=np.max(ctrsw),
                       extend='both')
    plt.contourf(x, y, sig1w, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax6.set_title('f', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax7 = fig.add_subplot(4, 2, 7)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot7 = m.contourf(x, y, winterr, ctrsr,
                       cmap=cmp1, vmin=np.min(ctrsr), vmax=np.max(ctrsr),
                       extend='both')
    plt.contourf(x, y, sig0r, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax7.set_title('g', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)
    ax8 = fig.add_subplot(4, 2, 8)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot8 = m.contourf(x, y, summerr, ctrsr,
                       cmap=cmp1, vmin=np.min(ctrsr), vmax=np.max(ctrsr),
                       extend='both')
    plt.contourf(x, y, sig1r, levels=[0, .5, 1.5],
                 hatches=["", "//"], alpha=0)
    ax8.set_title('h', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)

    cbar_ax = fig.add_axes([0.925, 0.7625, 0.005, 0.175])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference ($^{\circ}$C)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.535, 0.005, 0.175])
    b = fig.colorbar(plot4, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference ($^{\circ}$C)', size=20,
                fontsize=20, fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.3075, 0.005, 0.175])
    b = fig.colorbar(plot6, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference ($^{\circ}$C)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.925, 0.075, 0.005, 0.175])
    b = fig.colorbar(plot8, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$\mathregular{^{-1}}$)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 15  # in points
    ax1.annotate('DJF', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('JJA', xy=(0.5, 1), xytext=(0, pad), fontsize=32,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax1.annotate('T$\mathregular{_{mean}}$', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('TX', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax5.annotate('WBGT', xy=(0, 0.5), xytext=(-ax5.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax5.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax7.annotate('Precip', xy=(0, 0.5), xytext=(-ax7.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax7.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0.05, wspace=0, top=.95, bottom=0.05,
                        left=.05, right=.9)
