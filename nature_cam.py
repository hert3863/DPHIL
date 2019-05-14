# -*- coding: utf-8 -*-
"""
Created on Thu Aug  3 15:57:32 2017

@author: bakerh
"""


import numpy as np


def mainC():
    wbgt, wbgt95 = wbgt_calcC()
    tmeans = monthmeanTallC()
    tx = tmaxC()
    r95 = R95pC()
    r = dailymeanC('pr')
    tbar = dailymeanC('tas')
    #total_r, total_t, globprcnts = totalsC()
    #prcnts, all_mean = percentsC(tbar, tx, wbgt95, r, r95)
    tx90p = {}
    tx90p['All-Hist'] = tx['All-Hist'][2]
    tx90p['Plus15-Future'] = tx['Plus15-Future'][2]
    tx90p['Plus15-Future_LCO2'] = tx['Plus15-Future_LCO2'][2]
    a_t, b_t, c_t = mlrclusterC(tmeans, tx90p)
    a_w, b_w, c_w = mlrclusterC(tmeans, wbgt95)
    a_r, b_r, c_r = mlrclusterC(tmeans, r95, msk='TROP')
    return tbar, tmeans, tx, r, r95, wbgt95, a_t, b_t, c_t, a_w, b_w, c_w, a_r, b_r, c_r, prcnts, all_mean, total_r, total_t, globprcnts


def importallC():
    wbgt, wbgt95 = wbgt_calcC()
    tmeans = monthmeanTallC()
    tx = tmaxC()
    r95 = R95pC()
    r = dailymeanC('pr')
    tbar = dailymeanC('tas')
    return tbar, tmeans, tx, r, r95, wbgt95


def dailymeanC(item, region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}

    exps = ['All-Hist', 'Plus15-Future', 'Plus15-Future_LCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp + '/' +
                      item + '/*')
        rbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[920*e:920*(e+1)] = nc_fid.variables[item
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[900*e:900*(e+1)] = nc_fid.variables[item
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            rbar[0, y, :] = np.mean(wdata[900*y:900*(y+1), :], axis=0)
            rbar[1, y, :] = np.mean(sdata[920*y:920*(y+1), :], axis=0)

        output[exp] = rbar  # * 86400
    return output


def dailymean_meanC(item, region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    sindices = sindices.astype(int)
    output = {}

    exps = ['Plus15-Future-Future', 'Plus15-Future-Future_LCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree-2degree/' + exp + '/day/' +
                      item + '/*')
        sdata = np.zeros((region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata += np.mean(nc_fid.variables[item][sindices, 0, region[0]:
                                                    region[1],
                                                    region[2]:region[3]],
                             axis=0)
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))

        output[exp] = sdata * 86400 / len(a)
    return output


def monthlymeanC(item, level=9, region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*3))
    windices = np.zeros((10*3))
    for i in range(10):
        sindices[3*i:3*(i+1)] = np.linspace(5+12*i, 7+12*i, 3)
    for i in range(10):
        windices[3*i:3*(i+1)] = np.concatenate((np.linspace(12*i, 1+12*i,
                                               2), [11+12*i]))
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}

    exps = ['Plus15-Future', 'Plus15-Future_LCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp + '/mon/' +
                      item + '/*')
        rbar = np.zeros((2, np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        sdata = np.zeros((3*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((3*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[30*e:30*(e+1)] = nc_fid.variables[item
                                                    ][sindices, level, region[0]:
                                                      region[1],
                                                      region[2]:region[3]]
            wdata[30*e:30*(e+1)] = nc_fid.variables[item
                                                    ][windices, level, region[0]:
                                                      region[1],
                                                      region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            rbar[0, y, :] = np.mean(wdata[30*y:30*(y+1), :], axis=0)
            rbar[1, y, :] = np.mean(sdata[30*y:30*(y+1), :], axis=0)

        output[exp] = rbar  # * 86400
    return output


def R95pC(region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/All-Hist' +
                  '/day/pr/*')
    sdatar95p = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    wdatar95p = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatar95p[920*e:920*(e+1)] = nc_fid.variables['pr'
                                                      ][sindices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        wdatar95p[900*e:900*(e+1)] = nc_fid.variables['pr'
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
    exps = ['All-Hist', 'Plus15-Future', 'Plus15-Future_LCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/day/pr/*')

        sdata = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))

        r95p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                        region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[920*e:920*(e+1)] = nc_fid.variables['pr'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            wdata[900*e:900*(e+1)] = nc_fid.variables['pr'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]] * 86400
            print('Done precip ' + str(exp) + ' ' + str(e+1))

        truthw = np.zeros((90*10*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        truths = np.zeros((92*10*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(10*np.ma.size(a)):
            truthw[90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                        r95[0, :]).astype(int)
            truths[92*k:92*(k+1), :] = (sdata[92*k:92*(k+1), :] >
                                        r95[1, :]).astype(int)
        for m in range(np.ma.size(a)):
            r95p[0, m, :] = np.sum(truthw[900*m:900*(m+1), :], axis=0)
            r95p[1, m, :] = np.sum(truths[920*m:920*(m+1), :], axis=0)
        r95p[0, :] = r95p[0, :] / 10
        r95p[1, :] = r95p[1, :] / 10
        outputr95p[exp] = r95p  # percentage
    return outputr95p


def R95pC_un(region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/All-Hist' +
                  '/day/pr/*')
    sdatar95p = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    wdatar95p = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatar95p[920*e:920*(e+1)] = nc_fid.variables['pr'
                                                      ][sindices, 0,
                                                        region[0]: region[1],
                                                        region[2]:region[3]]
        wdatar95p[900*e:900*(e+1)] = nc_fid.variables['pr'
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
    exps = ['Plus15-Future', 'Plus15-Future_LCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/day/pr/*')

        sdata = np.zeros((92*10, region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*10, region[1]-region[0],
                         region[3]-region[2]))

        r95p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                        region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata = nc_fid.variables['pr'][sindices, 0, region[0]:
                                           region[1],
                                           region[2]:region[3]] * 86400
            wdata = nc_fid.variables['pr'][windices, 0, region[0]:
                                           region[1],
                                           region[2]:region[3]] * 86400
            print('Done precip ' + str(exp) + ' ' + str(e+1))

            truthw = np.zeros((90*10,
                              region[1]-region[0], region[3]-region[2]))
            truths = np.zeros((92*10,
                              region[1]-region[0], region[3]-region[2]))
            for k in range(10):
                truthw[90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                            r95[0, :]).astype(int)
                truths[92*k:92*(k+1), :] = (sdata[92*k:92*(k+1), :] >
                                            r95[1, :]).astype(int)

            r95p[0, e, :] = np.sum(truthw, axis=0)
            r95p[1, e, :] = np.sum(truths, axis=0)
        r95p[0, :] = r95p[0, :] / 10
        r95p[1, :] = r95p[1, :] / 10
        outputr95p[exp] = r95p  # percentage
    return outputr95p


def tmaxC(region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/All-Hist' +
                  '/day/tasmax/*')
    sdatatx90p = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    wdatatx90p = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatatx90p[920*e:920*(e+1)] = nc_fid.variables['tasmax'
                                                       ][sindices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        wdatatx90p[900*e:900*(e+1)] = nc_fid.variables['tasmax'
                                                       ][windices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        print('Done TX90p calc ' + str(e+1))
    daysw = np.arange(0, 900*np.ma.size(a), 90).astype(int)
    dayss = np.arange(0, 920*np.ma.size(a), 92).astype(int)
    t90w = np.zeros((90, region[1]-region[0],
                    region[3]-region[2]))
    t90s = np.zeros((92, region[1]-region[0],
                    region[3]-region[2]))
    for i in range(90):
        t90w[i, :] = np.percentile(wdatatx90p[daysw+i, :], 90, axis=0)
        t90s[i, :] = np.percentile(sdatatx90p[dayss+i, :], 90, axis=0)

    exps = ['All-Hist', 'Plus15-Future', 'Plus15-Future_LCO2']
    output = {}

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/day/tasmax/*')
        tx90p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                          region[3]-region[2]))
        txx = np.zeros((2, np.ma.size(a), 11, region[1]-region[0],
                        region[3]-region[2]))
        tx = np.zeros((2, np.ma.size(a), region[1]-region[0],
                       region[3]-region[2]))
        sdata = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[920*e:920*(e+1)] = nc_fid.variables['tasmax'
                                                      ][sindices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            wdata[900*e:900*(e+1)] = nc_fid.variables['tasmax'
                                                      ][windices, 0, region[0]:
                                                        region[1],
                                                        region[2]:region[3]]
            print('Done txmax ' + str(exp) + ' ' + str(e+1))

        truthw = np.zeros((90*10*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        truths = np.zeros((92*10*np.ma.size(a),
                          region[1]-region[0], region[3]-region[2]))
        for k in range(10*np.ma.size(a)):
            truthw[90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                        t90w).astype(int)
            truths[92*k:92*(k+1), :] = (sdata[92*k:92*(k+1), :] >
                                        t90s).astype(int)
        for m in range(np.ma.size(a)):
            tx90p[0, m, :] = np.sum(truthw[900*m:900*(m+1), :], axis=0)
            tx90p[1, m, :] = np.sum(truths[920*m:920*(m+1), :], axis=0)
        tx90p[0, :] = tx90p[0, :] / 10
        tx90p[1, :] = tx90p[1, :] / 10
        for y in range(np.ma.size(a)):
            tx[0, y, :] = np.mean(wdata[900*y:900*(y+1), :], axis=0)
            tx[1, y, :] = np.mean(sdata[920*y:920*(y+1), :], axis=0)
            for t in range(10):
                txx[0, y, t, :] = np.max(wdata[900*y+90*t:900*y+90*(t+1), :],
                                         axis=0)
                txx[1, y, t, :] = np.max(sdata[920*y+92*t:920*y+92*(t+1), :],
                                         axis=0)
        txx = np.mean(txx, axis=2)
        output[exp] = tx, txx, tx90p
    return output


def tmaxC_un(region=[0, 96, 0, 144]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/All-Hist' +
                  '/day/tasmax/*')
    sdatatx90p = np.zeros((92*10*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    wdatatx90p = np.zeros((90*10*np.ma.size(a), region[1]-region[0],
                           region[3]-region[2]))
    daysw = np.arange(0, 900*np.ma.size(a), 90).astype(int)
    dayss = np.arange(0, 920*np.ma.size(a), 92).astype(int)
    t90w = np.zeros((90, region[1]-region[0],
                    region[3]-region[2]))
    t90s = np.zeros((92, region[1]-region[0],
                    region[3]-region[2]))
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdatatx90p[920*e:920*(e+1)] = nc_fid.variables['tasmax'
                                                       ][sindices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        wdatatx90p[900*e:900*(e+1)] = nc_fid.variables['tasmax'
                                                       ][windices, 0,
                                                         region[0]: region[1],
                                                         region[2]:region[3]]
        print('Done TX90p calc ' + str(e+1))

    for i in range(90):
        t90w[i, :] = np.percentile(wdatatx90p[daysw+i, :], 90, axis=0)
    for i in range(92):
        t90s[i, :] = np.percentile(sdatatx90p[dayss+i, :], 90, axis=0)

    exps = ['Plus15-Future', 'Plus15-Future_LCO2']
    output = {}

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/day/tasmax/*')
        tx90p = np.zeros((2, np.ma.size(a), region[1]-region[0],
                          region[3]-region[2]))
        sdata = np.zeros((92*10, region[1]-region[0],
                         region[3]-region[2]))
        wdata = np.zeros((90*10, region[1]-region[0],
                         region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata = nc_fid.variables['tasmax'][sindices, 0, region[0]:
                                               region[1],
                                               region[2]:region[3]]
            wdata = nc_fid.variables['tasmax'][windices, 0, region[0]:
                                               region[1],
                                               region[2]:region[3]]
            print('Done txmax ' + str(exp) + ' ' + str(e+1))

            truthw = np.zeros((90*10,
                              region[1]-region[0], region[3]-region[2]))
            truths = np.zeros((92*10,
                              region[1]-region[0], region[3]-region[2]))
            for k in range(10):
                truthw[90*k:90*(k+1), :] = (wdata[90*k:90*(k+1), :] >
                                            t90w).astype(int)
                truths[92*k:92*(k+1), :] = (sdata[92*k:92*(k+1), :] >
                                            t90s).astype(int)
            tx90p[0, e, :] = np.sum(truthw, axis=0)
            tx90p[1, e, :] = np.sum(truths, axis=0)
        tx90p[0, :] = tx90p[0, :] / 10
        tx90p[1, :] = tx90p[1, :] / 10

        output[exp] = tx90p
    return output


def wbgt_calcC():
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    output = {}
    output_95p = {}

    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + 'All-Hist' +
                  '/day/tas/*')
    b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + 'All-Hist' +
                  '/day/hurs/*')
    wbgt = np.zeros((2, np.ma.size(a), 96, 144))
    wbgt_95p = np.zeros((2, np.ma.size(a), 96, 144))
    sdataT = np.zeros((3*10*np.ma.size(a), 96, 144))
    wdataT = np.zeros((3*10*np.ma.size(a), 96, 144))
    sdataH = np.zeros((3*10*np.ma.size(a), 96, 144))
    wdataH = np.zeros((3*10*np.ma.size(a), 96, 144))
    wbgts = np.zeros((3*10*np.ma.size(a), 96, 144))
    wbgtw = np.zeros((3*10*np.ma.size(a), 96, 144))
    jan = np.zeros((31*10))
    feb = np.zeros((28*10))
    jun = np.zeros((30*10))
    jul = np.zeros((31*10))
    aug = np.zeros((31*10))
    dec = np.zeros((31*10))
    for i in range(10):
        jan[31*i:31*(i+1)] = np.arange(0, 31) + i*90
        feb[28*i:28*(i+1)] = np.arange(0, 28) + i*90 + 31
        jun[30*i:30*(i+1)] = np.arange(0, 30) + i*92
        jul[31*i:31*(i+1)] = np.arange(0, 31) + i*92 + 30
        aug[31*i:31*(i+1)] = np.arange(0, 31) + i*92 + 61
        dec[31*i:31*(i+1)] = np.arange(0, 31) + i*90 + 59
    jan = jan.astype(int)
    feb = feb.astype(int)
    jun = jun.astype(int)
    jul = jul.astype(int)
    aug = aug.astype(int)
    dec = dec.astype(int)
    for e, d in enumerate(a):
        nc_fid = Dataset(d, 'r')
        sdataT_d = nc_fid.variables['tas'][sindices, 0, :]
        wdataT_d = nc_fid.variables['tas'][windices, 0, :]
        sdataT[30*e:30*(e+1):3] = np.mean(sdataT_d[jun, :], axis=0)
        sdataT[30*e+1:30*(e+1)+1:3] = np.mean(sdataT_d[jul, :], axis=0)
        sdataT[30*e+2:30*(e+1)+2:3] = np.mean(sdataT_d[aug, :], axis=0)
        wdataT[30*e:30*(e+1):3] = np.mean(wdataT_d[jan, :], axis=0)
        wdataT[30*e+1:30*(e+1)+1:3] = np.mean(wdataT_d[feb, :], axis=0)
        wdataT[30*e+2:30*(e+1)+2:3] = np.mean(wdataT_d[dec, :], axis=0)
        print('Done All-Hist tas ' + str(e+1))
    for e, d in enumerate(b):
        nc_fid = Dataset(d, 'r')
        sdataH_d = nc_fid.variables['hurs'][sindices, 0, :]
        wdataH_d = nc_fid.variables['hurs'][windices, 0, :]
        sdataH[30*e:30*(e+1):3] = np.mean(sdataH_d[jun, :], axis=0)
        sdataH[30*e+1:30*(e+1)+1:3] = np.mean(sdataH_d[jul, :], axis=0)
        sdataH[30*e+2:30*(e+1)+2:3] = np.mean(sdataH_d[aug, :], axis=0)
        wdataH[30*e:30*(e+1):3] = np.mean(wdataH_d[jan, :], axis=0)
        wdataH[30*e+1:30*(e+1)+1:3] = np.mean(wdataH_d[feb, :], axis=0)
        wdataH[30*e+2:30*(e+1)+2:3] = np.mean(wdataH_d[dec, :], axis=0)
        print('Done All-Hist hurs ' + str(e+1))

    Es = 6.112*np.exp(17.62*(sdataT-273.15)/(243.12+sdataT-273.15))
    Ew = 6.112*np.exp(17.62*(wdataT-273.15)/(243.12+wdataT-273.15))
    es = sdataH*Es/100
    ew = wdataH*Ew/100
    wbgts = 0.567*(sdataT-273.15)+0.393*es+3.94
    wbgtw = 0.567*(wdataT-273.15)+0.393*ew+3.94

    for y in range(np.ma.size(a)):
        wbgt[0, y, :] = np.mean(wbgtw[30*y:30*(y+1), :], axis=0)
        wbgt[1, y, :] = np.mean(wbgts[30*y:30*(y+1), :], axis=0)

    days = np.arange(0, 30*np.ma.size(a), 3).astype(int)
    months = np.arange(0, 30, 3).astype(int)
    wbgt95 = np.zeros((2, 30, 96, 144))
    for i in range(3):
        wbgt95[0, months+i, :] = np.percentile(wbgtw[days+i, :], 95, axis=0)
        wbgt95[1, months+i, :] = np.percentile(wbgts[days+i, :], 95, axis=0)

    for m in range(np.ma.size(a)):
        wbgt_95p[0, m] = np.sum((wbgtw[30*m:30*(m+1), :] - wbgt95[0]) *
                                (wbgtw[30*m:30*(m+1), :] >
                                wbgt95[0]).astype(int), axis=0)
        wbgt_95p[1, m] = np.sum((wbgts[30*m:30*(m+1), :] - wbgt95[1]) *
                                (wbgts[30*m:30*(m+1), :] >
                                wbgt95[1]).astype(int), axis=0)

    output['All-Hist'] = wbgt
    output_95p['All-Hist'] = wbgt_95p/10

    exps = ['Plus15-Future', 'Plus15-Future_LCO2']

    for x, ep in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + ep +
                      '/day/tas/*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + ep +
                      '/day/hurs/*')
        wbgt = np.zeros((2, np.ma.size(a), 96, 144))
        wbgt_95p = np.zeros((2, np.ma.size(a), 96, 144))
        sdataT = np.zeros((3*10*np.ma.size(a), 96, 144))
        wdataT = np.zeros((3*10*np.ma.size(a), 96, 144))
        sdataH = np.zeros((3*10*np.ma.size(a), 96, 144))
        wdataH = np.zeros((3*10*np.ma.size(a), 96, 144))
        wbgts = np.zeros((3*10*np.ma.size(a), 96, 144))
        wbgtw = np.zeros((3*10*np.ma.size(a), 96, 144))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdataT_d = nc_fid.variables['tas'][sindices, 0, :]
            wdataT_d = nc_fid.variables['tas'][windices, 0, :]
            sdataT[30*e:30*(e+1):3] = np.mean(sdataT_d[jun, :], axis=0)
            sdataT[30*e+1:30*(e+1)+1:3] = np.mean(sdataT_d[jul, :], axis=0)
            sdataT[30*e+2:30*(e+1)+2:3] = np.mean(sdataT_d[aug, :], axis=0)
            wdataT[30*e:30*(e+1):3] = np.mean(wdataT_d[jan, :], axis=0)
            wdataT[30*e+1:30*(e+1)+1:3] = np.mean(wdataT_d[feb, :], axis=0)
            wdataT[30*e+2:30*(e+1)+2:3] = np.mean(wdataT_d[dec, :], axis=0)
            print('Done tas' + str(ep) + ' ' + str(e+1))
        for e, d in enumerate(b):
            nc_fid = Dataset(d, 'r')
            sdataH_d = nc_fid.variables['hurs'][sindices, 0, :]
            wdataH_d = nc_fid.variables['hurs'][windices, 0, :]
            sdataH[30*e:30*(e+1):3] = np.mean(sdataH_d[jun, :], axis=0)
            sdataH[30*e+1:30*(e+1)+1:3] = np.mean(sdataH_d[jul, :], axis=0)
            sdataH[30*e+2:30*(e+1)+2:3] = np.mean(sdataH_d[aug, :], axis=0)
            wdataH[30*e:30*(e+1):3] = np.mean(wdataH_d[jan, :], axis=0)
            wdataH[30*e+1:30*(e+1)+1:3] = np.mean(wdataH_d[feb, :], axis=0)
            wdataH[30*e+2:30*(e+1)+2:3] = np.mean(wdataH_d[dec, :], axis=0)
            print('Done hurs ' + str(ep) + ' ' + str(e+1))

        Es = 6.112*np.exp(17.62*(sdataT-273.15)/(243.12+sdataT-273.15))
        Ew = 6.112*np.exp(17.62*(wdataT-273.15)/(243.12+wdataT-273.15))
        es = sdataH*Es/100
        ew = wdataH*Ew/100
        wbgts = 0.567*(sdataT-273.15)+0.393*es+3.94
        wbgtw = 0.567*(wdataT-273.15)+0.393*ew+3.94

        for y in range(np.ma.size(a)):
            wbgt[0, y, :] = np.mean(wbgtw[30*y:30*(y+1), :], axis=0)
            wbgt[1, y, :] = np.mean(wbgts[30*y:30*(y+1), :], axis=0)

        for m in range(np.ma.size(a)):
            wbgt_95p[0, m] = np.sum((wbgtw[30*m:30*(m+1), :] - wbgt95[0]) *
                                    (wbgtw[30*m:30*(m+1), :] >
                                    wbgt95[0]).astype(int), axis=0)
            wbgt_95p[1, m] = np.sum((wbgts[30*m:30*(m+1), :] - wbgt95[1]) *
                                    (wbgts[30*m:30*(m+1), :] >
                                    wbgt95[1]).astype(int), axis=0)
        output[ep] = wbgt
        output_95p[ep] = wbgt_95p/10
    return output, output_95p


def latweightmeanC(data, msk='no'):
    from netcdfread import ncread
    lsm = ncread('/home/bakerh/Documents/DPhil/Python/lsm_n72.nc', 'lsm')
    lsm1 = np.ones((96, 144))
    for i in range(96):
        for j in range(144):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    lat = np.linspace(-90, 90, 96)
    lon = np.linspace(0, 357.5, 144)
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    weighted = data * meshlatweight
    if msk == 'yes':
        mask = np.zeros((np.ma.size(data, axis=0), 96, 144))
        mask[:] = lsm1
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'NH':
        mask = np.zeros((np.ma.size(data, axis=0), 48, 144))
        lsm1 = lsm1[48:, :]
        mask[:] = lsm1
        weighted = weighted[:, 48:, :]
        meshlatweight = meshlatweight[48:, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'TROP':
        mask = np.zeros((np.ma.size(data, axis=0), 32, 144))
        lsm1 = lsm1[32:64, :]
        mask[:] = lsm1
        weighted = weighted[:, 32:64, :]
        meshlatweight = meshlatweight[32:64, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'MID':
        mask = np.zeros((np.ma.size(data, axis=0), 32, 144))
        lsm1 = lsm1[64:, :]
        mask[:] = lsm1
        weighted = weighted[:,  64:, :]
        meshlatweight = meshlatweight[64:, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'no':
        meaned = np.mean(weighted) / np.mean(meshlatweight)
    else:
        mask = np.zeros((np.ma.size(data, axis=0), msk[1]-msk[0],
                         msk[3]-msk[2]))
        lsm1 = lsm1[msk[0]:msk[1], msk[2]:msk[3]]
        mask[:] = lsm1
        weighted = weighted[:, msk[0]:msk[1], msk[2]:msk[3]]
        meshlatweight = meshlatweight[msk[0]:msk[1], msk[2]:msk[3]]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    return meaned


def totalsC():
    import glob
    from netCDF4 import Dataset
    exps = ['All-Hist', 'Plus15-Future', 'Plus15-Future_LCO2']
    total_r = np.zeros((3, 2, 3))  # exp; global, land only;
    total_t = np.zeros((3, 2, 3))  # all, winter, summer
    percentiles = np.zeros((4, 3, 5))
    globprcnts = np.zeros((4, 3, 5))
    '''
    sindices = np.zeros((10*92))
    windices = np.zeros((10*90))
    for i in range(10):
        sindices[92*i:92*(i+1)] = np.linspace(151+365*i, 242+365*i, 92)
    for i in range(10):
        windices[90*i:90*(i+1)] = np.concatenate((np.linspace(365*i, 58+365*i,
                                                 59), np.linspace(334+365*i,
                                                 364+365*i, 31)))
    sindices = sindices.astype(int)
    windices = windices.astype(int)
    '''
    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/pr/*')
        b = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/tas/*')
        dataR = np.zeros((365*10, 96, 144))
        dataT = np.zeros((365*10, 96, 144))
        dataR_11 = np.zeros((2, np.ma.size(a)))
        dataT_11 = np.zeros((2, np.ma.size(b)))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            dataR = nc_fid.variables['pr'][:, 0, :]
            dataR_11[0, e] = latweightmeanC(dataR, msk='no') * 86400
            dataR_11[1, e] = latweightmeanC(dataR, msk='yes') * 86400
            print('Done R totals ' + str(exp) + ' ' + str(e+1))
        #total_r[x, 0, 0] = latweightmeanC(dataR, msk='no') * 86400
        #total_r[x, 0, 1] = latweightmeanC(dataR[windices, :], msk='no') * 86400
        #total_r[x, 0, 2] = latweightmeanC(dataR[sindices, :], msk='no') * 86400
        #total_r[x, 1, 0] = latweightmeanC(dataR, msk='yes') * 86400
        #dataRw = dataR[windices, :]
        #dataRs = dataR[sindices, :]
        #total_r[x, 1, 1] = latweightmeanC(dataRw, msk='yes') * 86400
        #total_r[x, 1, 2] = latweightmeanC(dataRs, msk='yes') * 86400
        '''
        dataR_point = np.zeros((np.ma.size(a)*3960))
        for l in range(np.ma.size(a)*3960):
            dataR_point[l] = latweightmean(dataR[l, :]) * 86400
        dataR_11 = np.zeros((np.ma.size(a)))
        for y in range(np.ma.size(a)):
            dataR_11[y] = np.mean(dataR_point[3960*y:3960*(y+1)], axis=0)
        '''

        total_r[x, 0, 0] = np.mean(dataR_11[0, :])
        total_r[x, 1, 0] = np.mean(dataR_11[1, :])
        percentiles[2, x, :] = np.percentile(dataR_11[0, :], [10, 50, 90, 5, 95])
        percentiles[3, x, :] = np.percentile(dataR_11[1, :], [10, 50, 90, 5, 95])

        for f, g in enumerate(b):
            nc_fid = Dataset(g, 'r')
            dataT = nc_fid.variables['tas'][:, 0, :]
            dataT_11[0, f] = latweightmeanC(dataT, msk='no')
            dataT_11[1, f] = latweightmeanC(dataT, msk='yes')
            print('Done T totals ' + str(exp) + ' ' + str(f+1))
        #total_t[x, 0, 0] = latweightmeanC(dataT, msk='no')
        #total_t[x, 0, 1] = latweightmeanC(dataT[windices, :], msk='no')
        #total_t[x, 0, 2] = latweightmeanC(dataT[sindices, :], msk='no')
        #total_t[x, 1, 0] = latweightmeanC(dataT, msk='yes')
        #dataTw = dataT[windices, :]
        #dataTs = dataT[sindices, :]
        #total_t[x, 1, 1] = latweightmeanC(dataTw, msk='yes')
        #total_t[x, 1, 2] = latweightmeanC(dataTs, msk='yes')
        '''
        dataT_point = np.zeros((np.ma.size(b)*132))
        for l in range(np.ma.size(b)*132):
            dataT_point[l] = latweightmean(dataT[l, :])
        dataT_11 = np.zeros((np.ma.size(b)))
        for y in range(np.ma.size(b)):
            dataT_11[y] = np.mean(dataT_point[132*y:132*(y+1)], axis=0)
        '''

        total_t[x, 0, 0] = np.mean(dataT_11[0, :])
        total_t[x, 1, 0] = np.mean(dataT_11[1, :])
        percentiles[0, x, :] = np.percentile(dataT_11[0, :], [10, 50, 90, 5, 95])
        percentiles[1, x, :] = np.percentile(dataT_11[1, :], [10, 50, 90, 5, 95])
    for b in range(4):
        for a in range(3):
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


def monthmeanTallC(item):
    import glob
    from netCDF4 import Dataset

    output = {}

    exps = ['All-Hist', 'Plus15-Future', 'Plus15-Future_LCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/' + exp +
                      '/day/'+item+'/*')
        data = np.zeros((120, 96, 144))
        tbar = np.zeros((np.ma.size(a)))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            data = nc_fid.variables[item][:]
            print('Done tas ' + str(exp) + ' ' + str(e+1))
            tbar[e] = latweightmeanC(data, msk='no')
        output[exp] = tbar
    return output


def percentsC(tbar, txmax, wbgt_95p, r, r95p):
    exps = ['All-Hist', 'Plus15-Future', 'Plus15-Future_LCO2']
    percentiles = np.zeros((3, 10, 3))
    prcnts = np.zeros((3, 10, 3))
    all_mean = np.zeros((3, 10))

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
                means.append(latweightmeanC(np.expand_dims(var[j, :],
                                                          axis=0), msk='MID'))
                means_trop.append(latweightmeanC(np.expand_dims(var[j, :],
                                                               axis=0),
                                                msk='TROP'))
            percentiles[e, 2*i, :] = np.percentile(means, [10, 50, 90])
            percentiles[e, 2*i+1, :] = np.percentile(means_trop, [10, 50, 90])
            all_mean[e, 2*i] = np.mean(means)
            all_mean[e, 2*i+1] = np.mean(means_trop)

    for a in range(3):
        for b in range(10):
            prcnts[a, b, 0] = percentiles[a, b, 1]
            prcnts[a, b, 1] = percentiles[a, b, 1] - percentiles[a, b, 0]
            prcnts[a, b, 2] = percentiles[a, b, 2] - percentiles[a, b, 1]
    return prcnts, all_mean


def mlrclusterC(tmeans, varmeans, msk='MID'):
    import pandas as pd
    from statsmodels.formula.api import ols

    def elevenyrmns(data, msk='MID'):
        exps = ['Plus15-Future', 'Plus15-Future_LCO2']
        pmean = []
        for j in range(np.ma.size(data['All-Hist'][1], axis=0)):
            pmean.append(latweightmeanC(np.expand_dims(data['All-Hist']
                                                           [1, j, :],
                                                       axis=0), msk))
        pmean = np.mean(pmean)
        means = []
        size = []
        for exp in exps:
            size.append(np.ma.size(data[exp][1], axis=0))
            for j in range(np.ma.size(data[exp][1], axis=0)):
                means.append(latweightmeanC(np.expand_dims(data[exp][1, j, :],
                                                           axis=0), msk)-pmean)
        return means, size
    co2 = np.array([390.4, 423.4, 379.0])
    f = 3.71 * np.log(co2/278) / np.log(2)  # 278 as preindust [CO2]
    d_f = f[1:] - f[0]
    t_means = (np.concatenate((tmeans['Plus15-Future'], tmeans['Plus15-Future_LCO2'])) -
               np.mean(tmeans['All-Hist']))
    t_size = np.ma.size(t_means)
    var_means, v_size = elevenyrmns(varmeans, msk)
    if t_size == np.sum(v_size):
        d_fs = np.zeros((np.sum(v_size)))
        d_fs[:v_size[0]] = d_f[0]
        d_fs[v_size[0]:] = d_f[1]

    else:
        print("Error, Tmean and extreme differ in size")
    data = pd.DataFrame({'d_f': d_fs, 'd_t': t_means, 'd_var': var_means})
    model = ols("d_var ~ d_t + d_f - 1", data).fit()
    # print(model.summary())
    # print(model._results.params)
    beta, alpha = model._results.params
    cvs = model.cov_params()
    return alpha, beta, cvs


def mlrcluster_gridC(tmeans, varmeans):
    import pandas as pd
    from statsmodels.formula.api import ols
    beta = np.zeros((2, 96, 144))
    co2 = np.array([390.4, 423.4, 379.0])
    f = 3.71 * np.log(co2/278) / np.log(2)  # 278 as preindust [CO2]
    d_f = f[1:] - f[0]
    t_means = (np.concatenate((tmeans['Plus15-Future'],
                               tmeans['Plus15-Future_LCO2'])) -
               np.mean(tmeans['All-Hist']))
    t_size = np.ma.size(t_means)
    d_fs = np.zeros((np.sum(t_size)))
    d_fs[:np.ma.size(varmeans['Plus15-Future'][1], axis=0)] = d_f[0]
    d_fs[np.ma.size(varmeans['Plus15-Future'][1], axis=0):] = d_f[1]
    for i in range(96):
        for j in range(144):
            for k in range(2):
                var = np.concatenate((varmeans['Plus15-Future'][k, :, i, j],
                                      varmeans['Plus15-Future_LCO2'][k, :, i, j])) - np.mean(varmeans['All-Hist'][k, :, i, j], axis=0)
                data = pd.DataFrame({'d_f': d_fs, 'd_t': t_means,
                                     'd_var': var})
                model = ols("d_var ~ d_t + d_f - 1", data).fit()
                beta[k, i, j] = model._results.params[0]
        print(str(i+1) + '/96')
    return beta


def fig4C(tmeans, conversion, rcpplume, alpha, beta, cvs):
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
    co2 = np.array([390.4, 423.4, 379.0])
    tmean = np.zeros((3))
    tmean[0] = np.mean(tmeans['All-Hist'])
    tmean[1] = np.mean(tmeans['Plus15-Future'])
    tmean[2] = np.mean(tmeans['Plus15-Future_LCO2'])

    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv['d_f'][1] +
                       (x*alpha/beta**2)**2*cv['d_t'][0] -
                       2*(x**2*alpha*cv['d_t'][1]**2/beta**3))
        return sigT

    co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301

    ax1.set_ylim([0, 1.5])
    hfont = {'fontname': 'Arial'}
    ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2015 \
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
    plt.title('CAM4-2degree WBGT95p', y=1.01, fontsize=20, **hfont)
    # plt.subplots_adjust(hspace=0, wspace=0.05, top=.95, bottom=0.1, left=.05,
    #                     right=.97)
    un = oub - lb
    p_ch = np.array([float(ub)-oub, float(ub)-float(ubl),
                     float(ubu)-float(ub)])*100/un

    return oub, float(ub), float(ubl), float(ubu), p_ch

'''
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
    globprcnts3[1, 1:, 0] -= globprcnts[1, 0, 0]
    globprcnts3 = np.vstack((globprcnts3[1, 2], globprcnts3[1, 1],
                             globprcnts3[1, 3]))
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
'''


def fig1_daysC(datat, dataw, datar, betat=np.zeros(2), betaw=np.zeros(2),
               betar=np.zeros(2), rem='no'):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lat = np.linspace(-90, 90, 96)
    lon = np.linspace(0, 357.5, 144)
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    winterr = np.mean(datar['Plus15-Future'][0, :], axis=0) - np.mean(datar['Plus15-Future_LCO2'][0, :], axis=0)
    summerr = np.mean(datar['Plus15-Future'][1, :], axis=0) - np.mean(datar['Plus15-Future_LCO2'][1, :], axis=0)

    wintert = np.mean(datat['Plus15-Future'][0, :], axis=0) - np.mean(datat['Plus15-Future_LCO2'][0, :], axis=0)
    summert = np.mean(datat['Plus15-Future'][1, :], axis=0) - np.mean(datat['Plus15-Future_LCO2'][1, :], axis=0)

    winterw = np.mean(dataw['Plus15-Future'][0, :], axis=0) - np.mean(dataw['Plus15-Future_LCO2'][0, :], axis=0)
    summerw = np.mean(dataw['Plus15-Future'][1, :], axis=0) - np.mean(dataw['Plus15-Future_LCO2'][1, :], axis=0)

    dT = 0.038931112301099802
    if rem == 'yes':

        winterr -= dT * betar[0]
        summerr -= dT * betar[1]
        wintert -= dT * betat[0]
        summert -= dT * betat[1]
        winterw -= dT * betaw[0]
        summerw -= dT * betaw[1]

    sig0r = ttestC(datar['Plus15-Future_LCO2'][0, :],
                  datar['Plus15-Future'][0, :]-betar[0]*dT, siglevel=10)
    sig1r = ttestC(datar['Plus15-Future_LCO2'][1, :],
                  datar['Plus15-Future'][1, :]-betar[1]*dT, siglevel=10)
    sig0t = ttestC(datat['Plus15-Future'][0, :]-betat[0]*dT,
                  datat['Plus15-Future_LCO2'][0, :], siglevel=10)
    sig1t = ttestC(datat['Plus15-Future'][1, :]-betat[1]*dT,
                  datat['Plus15-Future_LCO2'][1, :], siglevel=10)
    sig0w = ttestC(dataw['Plus15-Future'][0, :]-betaw[0]*dT,
                  dataw['Plus15-Future_LCO2'][0, :], siglevel=10)
    sig1w = ttestC(dataw['Plus15-Future'][1, :]-betaw[1]*dT,
                  dataw['Plus15-Future_LCO2'][1, :], siglevel=10)


    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    lsm = ncread('/home/bakerh/Documents/DPhil/Python/lsm_n72.nc', 'lsm')
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
    ctrst = np.linspace(-7, 7, 17)
    ctrsw = np.linspace(-.4, .4, 17)
    ctrsr = np.linspace(-.8, .8, 17)
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
                 hatches=["", ".."], alpha=0, color='white')
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
                 hatches=["", ".."], alpha=0, color='white')
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
                 hatches=["", ".."], alpha=0, edgecolor='white')
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
                 hatches=["", ".."], alpha=0, edgecolor='white')
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
                 hatches=["", ".."], alpha=0, edgecolor='white')
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
                 hatches=["", ".."], alpha=0, color='white')
    ax6.set_title('f', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.9)

    cbar_ax = fig.add_axes([0.905, 0.7125, 0.005, 0.2])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (days season$\mathregular{^{-1}}$)', size=20, fontsize=20,
                fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.905, 0.405, 0.005, 0.2])
    b = fig.colorbar(plot4, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference ($^{\circ}$C months)', size=20,
                fontsize=20, fontname='Arial')
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=20)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.905, 0.1, 0.005, 0.2])
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


def ttestC(series1, series2, siglevel=10, testtype='two'):
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
    sigarray = np.full((96, 144), siglevel)
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


def kstestC(series1, series2, siglevel=10):
    import scipy.stats as st
    sigarray = np.full((96, 144), siglevel)
    p1 = np.zeros((96, 144))
    for i in range(96):
        for j in range(144):
            p1[i, j] = st.ks_2samp(series1[:, i, j], series2[:, i, j])[1]
    sig = np.greater_equal(sigarray, p1*100).astype(int)
    return sig


def fig1_extraC(datat, datatx, dataw, datar):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lat = np.linspace(-90, 90, 96)
    lon = np.linspace(0, 357.5, 144)
    #lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_522/atmos/item15201_monthly_mean/item15201_monthly_mean_a00a_2090-01_2100-12.nc', 'latitude1')
    winterr = np.mean(datar['Plus15-Future'][0, :], axis=0) - np.mean(datar['Plus15-Future_LCO2'][0, :], axis=0)
    summerr = np.mean(datar['Plus15-Future'][1, :], axis=0) - np.mean(datar['Plus15-Future_LCO2'][1, :], axis=0)
    sig0r = ttestC(datar['Plus15-Future_LCO2'][0, :],
                  datar['Plus15-Future'][0, :], siglevel=10)
    sig1r = ttestC(datar['Plus15-Future_LCO2'][1, :],
                  datar['Plus15-Future'][1, :], siglevel=10)

    wintert = np.mean(datat['Plus15-Future'][0, :], axis=0) - np.mean(datat['Plus15-Future_LCO2'][0, :], axis=0)
    summert = np.mean(datat['Plus15-Future'][1, :], axis=0) - np.mean(datat['Plus15-Future_LCO2'][1, :], axis=0)
    sig0t = ttestC(datat['Plus15-Future'][0, :],
                  datat['Plus15-Future_LCO2'][0, :], siglevel=10)
    sig1t = ttestC(datat['Plus15-Future'][1, :],
                  datat['Plus15-Future_LCO2'][1, :], siglevel=10)

    winterw = np.mean(dataw['Plus15-Future'][0, :], axis=0) - np.mean(dataw['Plus15-Future_LCO2'][0, :], axis=0)
    summerw = np.mean(dataw['Plus15-Future'][1, :], axis=0) - np.mean(dataw['Plus15-Future_LCO2'][1, :], axis=0)
    sig0w = ttestC(dataw['Plus15-Future'][0, :],
                  dataw['Plus15-Future_LCO2'][0, :], siglevel=10)
    sig1w = ttestC(dataw['Plus15-Future'][1, :],
                  dataw['Plus15-Future_LCO2'][1, :], siglevel=10)

    wintertx = np.mean(datatx['Plus15-Future'][0][0, :], axis=0) - np.mean(datatx['Plus15-Future_LCO2'][0][0, :], axis=0)
    summertx = np.mean(datatx['Plus15-Future'][0][1, :], axis=0) - np.mean(datatx['Plus15-Future_LCO2'][0][1, :], axis=0)
    sig0tx = ttestC(datatx['Plus15-Future'][0][0, :],
                  datatx['Plus15-Future_LCO2'][0][0, :], siglevel=10)
    sig1tx = ttestC(datatx['Plus15-Future'][0][1, :],
                  datatx['Plus15-Future_LCO2'][0][1, :], siglevel=10)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    lsm = ncread('/home/bakerh/Documents/DPhil/Python/lsm_n72.nc', 'lsm')
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
    ctrst = np.linspace(-.6, .6, 17)
    ctrsw = np.linspace(-.4, .4, 17)
    ctrsr = np.linspace(-.4, .4, 17)
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


def plotallC(data, colorlimit, meaning='mean', mask='yes', precip='no',
            cbarleg='cbarlg', pltlbl='pltlbl', b1='Plus15-Future_LCO2', b2='Plus15-Future'):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/All-Hist/va/va_Amon_CAM4-2degree-2degree_All-Hist_est1_v1-0_ens0000_200601-201512.nc', 'lon')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/CAM4-2degree/All-Hist/va/va_Amon_CAM4-2degree-2degree_All-Hist_est1_v1-0_ens0000_200601-201512.nc', 'lat')
    if meaning == 'mean':
        winter = np.mean(data[b2][0, :], axis=0) - np.mean(data[b1][0, :], axis=0)
        summer = np.mean(data[b2][1, :], axis=0) - np.mean(data[b1][1, :], axis=0)
        sig0 = ttestC(data[b1][0, :],
                     data[b2][0, :], siglevel=10)
        sig1 = ttestC(data[b1][1, :],
                     data[b2][1, :], siglevel=10)
    elif meaning == 'non':
        winter = data[b2][0, :] - data[b1][0, :]
        summer = data[b2][1, :] - data[b1][1, :]
        sig0 = np.zeros((145, 192))
        sig1 = np.zeros((145, 192))
    else:
        winter = np.mean(data[b2][meaning][0, :], axis=0) - np.mean(data[b1][meaning][0, :], axis=0)
        summer = np.mean(data[b2][meaning][1, :], axis=0) - np.mean(data[b1][meaning][1, :], axis=0)
        sig0 = ttestC(data[b2][meaning][0, :],
                     data[b1][meaning][0, :], siglevel=10)
        sig1 = ttestC(data[b2][meaning][1, :],
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
