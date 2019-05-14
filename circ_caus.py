#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 09:31:27 2017

@author: bakerh
"""


def ensids(vrs=['tasmax', 'pr', 'zg', 'rsds', 'rlds', 'rlds_rlus']):
    import glob
    import fnmatch
    ncs_H = {}
    ncs_L = {}
    ids_H = []
    ids_L = []
    for var in vrs:
        hco2 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/' + var + '/*')
        lco2 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/' + var + '/*')
        for i in range(len(hco2)):
            hco2[i] = hco2[i][-23:-19]
        for i in range(len(lco2)):
            lco2[i] = lco2[i][-23:-19]
        ncs_H[var] = hco2
        ncs_L[var] = lco2
    for item in (ncs_H[vrs[0]]):
        c = 0
        for j in range(len(vrs)-1):
            if fnmatch.filter(ncs_H[vrs[j+1]], item):
                c += 1
        if c == len(vrs)-1:
            ids_H.append(item)
    for item in (ncs_L[vrs[0]]):
        c = 0
        for j in range(len(vrs)-1):
            if fnmatch.filter(ncs_L[vrs[j+1]], item):
                c += 1
        if c == len(vrs)-1:
            ids_L.append(item)
    return ids_H, ids_L


def dailymean(item, code, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    import numpy as np
    sindices = np.zeros((11*90))
    for i in range(11):
        sindices[90*i:90*(i+1)] = np.linspace(150+360*i, 239+360*i, 90)
    sindices = sindices.astype(int)
    output = {}

    exps = ['Plus15-Future_LCO2', 'Plus15-Future_HCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                      exp + '/day/' + item + '/*')
        rbar = np.zeros((region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata = nc_fid.variables[code][sindices, 0, region[0]: region[1],
                                           region[2]:region[3]]

            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
            rbar += np.mean(sdata, axis=0)

        output[exp] = rbar / len(a)  # * 86400
    return output


def monthlymean(item, code, lev, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    import numpy as np
    sindices = np.zeros((11*3))
    for i in range(11):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    sindices = sindices.astype(int)
    output = {}

    exps = ['Plus15-Future_LCO2', 'Plus15-Future_HCO2']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                      exp + '/mon/' + item + '/*')
        rbar = np.zeros((33, region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata = nc_fid.variables[code][sindices, lev, region[0]: region[1],
                                           region[2]:region[3]]

            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
            rbar += sdata

        output[exp] = rbar / len(a)  # * 86400
    return output


def dailymeanSU(item, code, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    import numpy as np
    output = {}

    exps = ['Plus15-Future_LCO2_SU', 'Plus15-Future_HCO2_SU']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                      exp + '/day/' + item + '/*')
        rbar = np.zeros((1110, region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata = nc_fid.variables[code][:, 0, region[0]: region[1],
                                           region[2]:region[3]]

            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
            rbar += sdata

        output[exp] = rbar / len(a)  # * 86400
    return output


def monthlymeanSU(item, code, level, region=[0, 145, 0, 192]):
    import glob
    from netCDF4 import Dataset
    import numpy as np
    output = {}

    exps = ['Plus15-Future_LCO2_SU', 'Plus15-Future_HCO2_SU']

    for x, exp in enumerate(exps):
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                      exp + '/mon/' + item + '/*')
        rbar = np.zeros((37, region[1]-region[0],
                         region[3]-region[2]))

        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata = nc_fid.variables[code][:, level, region[0]: region[1],
                                           region[2]:region[3]]

            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
            rbar += sdata

        output[exp] = rbar / len(a)  # * 86400
    return output


def imports(ids_H, ids_L):
    import numpy as np
    from netcdfread import ncread
    tasmax_rus_H = np.zeros((len(ids_H), 3960))
    tasmax_eur_H = np.zeros((len(ids_H), 3960))
    tasmax_nus_H = np.zeros((len(ids_H), 3960))
    rsds_rus_H = np.zeros((len(ids_H), 3960))
    rsds_eur_H = np.zeros((len(ids_H), 3960))
    rsds_nus_H = np.zeros((len(ids_H), 3960))
    rlus_rus_H = np.zeros((len(ids_H), 3960))
    rlus_eur_H = np.zeros((len(ids_H), 3960))
    rlus_nus_H = np.zeros((len(ids_H), 3960))
    zg_rus_H = np.zeros((len(ids_H), 3960))
    zg_eur_H = np.zeros((len(ids_H), 3960))
    zg_nus_H = np.zeros((len(ids_H), 3960))
    pr_nsa_H = np.zeros((len(ids_H), 3960))
    pr_car_H = np.zeros((len(ids_H), 3960))
    pr_naf_H = np.zeros((len(ids_H), 3960))
    pr_scs_H = np.zeros((len(ids_H), 3960))
    pr_ind_H = np.zeros((len(ids_H), 3960))
    tasmax_rus_L = np.zeros((len(ids_L), 3960))
    tasmax_eur_L = np.zeros((len(ids_L), 3960))
    tasmax_nus_L = np.zeros((len(ids_L), 3960))
    rsds_rus_L = np.zeros((len(ids_L), 3960))
    rsds_eur_L = np.zeros((len(ids_L), 3960))
    rsds_nus_L = np.zeros((len(ids_L), 3960))
    rlus_rus_L = np.zeros((len(ids_L), 3960))
    rlus_eur_L = np.zeros((len(ids_L), 3960))
    rlus_nus_L = np.zeros((len(ids_L), 3960))
    zg_rus_L = np.zeros((len(ids_L), 3960))
    zg_eur_L = np.zeros((len(ids_L), 3960))
    zg_nus_L = np.zeros((len(ids_L), 3960))
    pr_nsa_L = np.zeros((len(ids_L), 3960))
    pr_car_L = np.zeros((len(ids_L), 3960))
    pr_naf_L = np.zeros((len(ids_L), 3960))
    pr_scs_L = np.zeros((len(ids_L), 3960))
    pr_ind_L = np.zeros((len(ids_L), 3960))

    for i in range(len(ids_L)):
        tasmax_rus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/tasmax/item3236_daily_maximum_' + ids_L[i] + '_2090-01_2100-12.nc', 'item3236_daily_maximum')[:, 0, 16:33, 37:65], axis=(1, 2))
        tasmax_eur_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/tasmax/item3236_daily_maximum_' + ids_L[i] + '_2090-01_2100-12.nc', 'item3236_daily_maximum')[:, 0, 28:41, 3:17], axis=(1, 2))
        tasmax_nus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/tasmax/item3236_daily_maximum_' + ids_L[i] + '_2090-01_2100-12.nc', 'item3236_daily_maximum')[:, 0, 24:37, 128:145], axis=(1, 2))
        rsds_rus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rsds/item1235_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item1235_daily_mean')[:, 0, 16:33, 37:65], axis=(1, 2))
        rsds_eur_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rsds/item1235_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item1235_daily_mean')[:, 0, 28:41, 3:17], axis=(1, 2))
        rsds_nus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rsds/item1235_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item1235_daily_mean')[:, 0, 24:37, 128:145], axis=(1, 2))
        rlus_rus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rlds/item2207_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item2207_daily_mean')[:, 0, 16:33, 37:65], axis=(1, 2)) - np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rlds_rlus/item2201_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item2201_daily_mean')[:, 0, 16:33, 37:65], axis=(1, 2))
        rlus_eur_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rlds/item2207_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item2207_daily_mean')[:, 0, 28:41, 3:17], axis=(1, 2)) - np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rlds_rlus/item2201_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item2201_daily_mean')[:, 0, 28:41, 3:17], axis=(1, 2))
        rlus_nus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rlds/item2207_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item2207_daily_mean')[:, 0, 24:37, 128:145], axis=(1, 2)) - np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/rlds_rlus/item2201_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item2201_daily_mean')[:, 0, 24:37, 128:145], axis=(1, 2))
        zg_rus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/zg/item16202_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item16202_daily_mean')[:, 0, 16:37, 48:70], axis=(1, 2))
        zg_eur_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/zg/item16202_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item16202_daily_mean')[:, 0, 20:37, 5:17], axis=(1, 2))
        zg_nus_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/zg/item16202_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item16202_daily_mean')[:, 0, 28:41, 157:171], axis=(1, 2))
        pr_nsa_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/pr/item5216_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 68:77, 149:161], axis=(1, 2))
        pr_car_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/pr/item5216_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 48:65, 139:166], axis=(1, 2))
        pr_naf_L[i] = np.nanmean(np.concatenate((ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/pr/item5216_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 60:65, 176:], ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/pr/item5216_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 60:65, :27]), axis=2), axis=(1, 2))
        pr_scs_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/pr/item5216_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 56:65, 59:75], axis=(1, 2))
        pr_ind_L[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/pr/item5216_daily_mean_' + ids_L[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 68:81, 53:75], axis=(1, 2))
        print('Done ids_l: ' + str(i))
    tasmax_rus_L = np.mean(np.mean(np.reshape(tasmax_rus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    tasmax_eur_L = np.mean(np.mean(np.reshape(tasmax_eur_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    tasmax_nus_L = np.mean(np.mean(np.reshape(tasmax_nus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    rsds_rus_L = np.mean(np.mean(np.reshape(rsds_rus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    rsds_eur_L = np.mean(np.mean(np.reshape(rsds_eur_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    rsds_nus_L = np.mean(np.mean(np.reshape(rsds_nus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    rlus_rus_L = np.mean(np.mean(np.reshape(rlus_rus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    rlus_eur_L = np.mean(np.mean(np.reshape(rlus_eur_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    rlus_nus_L = np.mean(np.mean(np.reshape(rlus_nus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    zg_rus_L = np.mean(np.mean(np.reshape(zg_rus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    zg_eur_L = np.mean(np.mean(np.reshape(zg_eur_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    zg_nus_L = np.mean(np.mean(np.reshape(zg_nus_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    pr_nsa_L = np.mean(np.mean(np.reshape(pr_nsa_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    pr_car_L = np.mean(np.mean(np.reshape(pr_car_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    pr_naf_L = np.mean(np.mean(np.reshape(pr_naf_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    pr_scs_L = np.mean(np.mean(np.reshape(pr_scs_L, (len(ids_L), 792, 5)), axis=2), axis=0)
    pr_ind_L = np.mean(np.mean(np.reshape(pr_ind_L, (len(ids_L), 792, 5)), axis=2), axis=0)

    for i in range(len(ids_H)):
        tasmax_rus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/tasmax/item3236_daily_maximum_' + ids_H[i] + '_2090-01_2100-12.nc', 'item3236_daily_maximum')[:, 0, 16:33, 37:65], axis=(1, 2))
        tasmax_eur_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/tasmax/item3236_daily_maximum_' + ids_H[i] + '_2090-01_2100-12.nc', 'item3236_daily_maximum')[:, 0, 28:41, 3:17], axis=(1, 2))
        tasmax_nus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/tasmax/item3236_daily_maximum_' + ids_H[i] + '_2090-01_2100-12.nc', 'item3236_daily_maximum')[:, 0, 24:37, 128:145], axis=(1, 2))
        rsds_rus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rsds/item1235_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item1235_daily_mean')[:, 0, 16:33, 37:65], axis=(1, 2))
        rsds_eur_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rsds/item1235_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item1235_daily_mean')[:, 0, 28:41, 3:17], axis=(1, 2))
        rsds_nus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rsds/item1235_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item1235_daily_mean')[:, 0, 24:37, 128:145], axis=(1, 2))
        rlus_rus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rlds/item2207_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item2207_daily_mean')[:, 0, 16:33, 37:65], axis=(1, 2)) - np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rlds_rlus/item2201_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item2201_daily_mean')[:, 0, 16:33, 37:65], axis=(1, 2))
        rlus_eur_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rlds/item2207_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item2207_daily_mean')[:, 0, 28:41, 3:17], axis=(1, 2)) - np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rlds_rlus/item2201_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item2201_daily_mean')[:, 0, 28:41, 3:17], axis=(1, 2))
        rlus_nus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rlds/item2207_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item2207_daily_mean')[:, 0, 24:37, 128:145], axis=(1, 2)) - np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/rlds_rlus/item2201_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item2201_daily_mean')[:, 0, 24:37, 128:145], axis=(1, 2))
        zg_rus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/zg/item16202_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item16202_daily_mean')[:, 0, 16:37, 48:70], axis=(1, 2))
        zg_eur_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/zg/item16202_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item16202_daily_mean')[:, 0, 20:37, 5:17], axis=(1, 2))
        zg_nus_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/zg/item16202_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item16202_daily_mean')[:, 0, 28:41, 157:171], axis=(1, 2))
        pr_nsa_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/pr/item5216_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 68:77, 149:161], axis=(1, 2))
        pr_car_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/pr/item5216_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 48:65, 139:166], axis=(1, 2))
        pr_naf_H[i] = np.nanmean(np.concatenate((ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/pr/item5216_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 60:65, 176:], ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/pr/item5216_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 60:65, :27]), axis=2), axis=(1, 2))
        pr_scs_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/pr/item5216_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 56:65, 59:75], axis=(1, 2))
        pr_ind_H[i] = np.nanmean(ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_HCO2/day/pr/item5216_daily_mean_' + ids_H[i] + '_2090-01_2100-12.nc', 'item5216_daily_mean')[:, 0, 68:81, 53:75], axis=(1, 2))
        print('Done ids_h: ' + str(i))
    tasmax_rus_H = np.mean(np.reshape(tasmax_rus_H, (len(ids_H), 792, 5)), axis=2)
    tasmax_eur_H = np.mean(np.reshape(tasmax_eur_H, (len(ids_H), 792, 5)), axis=2)
    tasmax_nus_H = np.mean(np.reshape(tasmax_nus_H, (len(ids_H), 792, 5)), axis=2)
    rsds_rus_H = np.mean(np.reshape(rsds_rus_H, (len(ids_H), 792, 5)), axis=2)
    rsds_eur_H = np.mean(np.reshape(rsds_eur_H, (len(ids_H), 792, 5)), axis=2)
    rsds_nus_H = np.mean(np.reshape(rsds_nus_H, (len(ids_H), 792, 5)), axis=2)
    rlus_rus_H = np.mean(np.reshape(rlus_rus_H, (len(ids_H), 792, 5)), axis=2)
    rlus_eur_H = np.mean(np.reshape(rlus_eur_H, (len(ids_H), 792, 5)), axis=2)
    rlus_nus_H = np.mean(np.reshape(rlus_nus_H, (len(ids_H), 792, 5)), axis=2)
    zg_rus_H = np.mean(np.reshape(zg_rus_H, (len(ids_H), 792, 5)), axis=2)
    zg_eur_H = np.mean(np.reshape(zg_eur_H, (len(ids_H), 792, 5)), axis=2)
    zg_nus_H = np.mean(np.reshape(zg_nus_H, (len(ids_H), 792, 5)), axis=2)
    pr_nsa_H = np.mean(np.reshape(pr_nsa_H, (len(ids_H), 792, 5)), axis=2)
    pr_car_H = np.mean(np.reshape(pr_car_H, (len(ids_H), 792, 5)), axis=2)
    pr_naf_H = np.mean(np.reshape(pr_naf_H, (len(ids_H), 792, 5)), axis=2)
    pr_scs_H = np.mean(np.reshape(pr_scs_H, (len(ids_H), 792, 5)), axis=2)
    pr_ind_H = np.mean(np.reshape(pr_ind_H, (len(ids_H), 792, 5)), axis=2)

    tasmax_rus = np.reshape(tasmax_rus_H - tasmax_rus_L,(-1))
    tasmax_eur = np.reshape(tasmax_eur_H - tasmax_eur_L,(-1))
    tasmax_nus = np.reshape(tasmax_nus_H - tasmax_nus_L,(-1))
    rsds_rus = np.reshape(rsds_rus_H - rsds_rus_L,(-1))
    rsds_eur = np.reshape(rsds_eur_H - rsds_eur_L,(-1))
    rsds_nus = np.reshape(rsds_nus_H - rsds_nus_L,(-1))
    rlus_rus = np.reshape(rlus_rus_H - rlus_rus_L,(-1))
    rlus_eur = np.reshape(rlus_eur_H - rlus_eur_L,(-1))
    rlus_nus = np.reshape(rlus_nus_H - rlus_nus_L,(-1))
    zg_rus = np.reshape(zg_rus_H - zg_rus_L,(-1))
    zg_eur = np.reshape(zg_eur_H - zg_eur_L,(-1))
    zg_nus = np.reshape(zg_nus_H - zg_nus_L,(-1))
    pr_nsa = np.reshape(pr_nsa_H - pr_nsa_L,(-1)) * 86400
    pr_car = np.reshape(pr_car_H - pr_car_L,(-1)) * 86400
    pr_naf = np.reshape(pr_naf_H - pr_naf_L,(-1)) * 86400
    pr_scs = np.reshape(pr_scs_H - pr_scs_L,(-1)) * 86400
    pr_ind = np.reshape(pr_ind_H - pr_ind_L,(-1)) * 86400

    data = np.transpose(np.vstack((tasmax_rus, tasmax_eur, tasmax_nus, rsds_rus, rsds_eur, rsds_nus, rlus_rus, rlus_eur, rlus_nus, zg_rus, zg_eur, zg_nus, pr_nsa, pr_car, pr_naf, pr_scs, pr_ind)))
    var_names = ['tasmax_rus', 'tasmax_eur', 'tasmax_nus', 'rsds_rus', 'rsds_eur', 'rsds_nus', 'rlus_rus', 'rlus_eur', 'rlus_nus', 'zg_rus', 'zg_eur', 'zg_nus', 'pr_nsa', 'pr_car', 'pr_naf', 'pr_scs', 'pr_ind']
    return data, var_names


def maplot(data, colormax=1, colormin=-999, mask='no', title='', precip='no'):
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
    pdata = data['Plus15-Future_HCO2'] - data['Plus15-Future_LCO2']
    #pdata = np.mean(data['Plus15-Future_HCO2_SU'][-210:-120], axis=0) - np.mean(data['Plus15-Future_LCO2_SU'][-210:-120], axis=0)
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
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
    b.set_label(label='U200 difference (hPa)')
    parallels = m.drawparallels(np.arange(-90., 91., 5.))
    meridians = m.drawmeridians(np.arange(-180., 181., 10))
    m.drawparallels(parallels, labels=[True, True, True, True], labelstyle='+/-')
    m.drawmeridians(meridians, labels=[True, True, True, True], labelstyle='+/-')
    plt.title(title, y=1.08)
    plt.subplots_adjust(hspace=0.05, wspace=0, top=.95, bottom=0.05,
                        left=.05, right=.95)
    plt.show()


def caus(data, var_names):
    import numpy as np
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    import sklearn

    import tigramite
    from tigramite import data_processing as pp
    from tigramite import plotting as tp
    from tigramite.pcmci import PCMCI
    from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb
    from tigramite.models import LinearMediation, Prediction

    data_mask_row = np.zeros(len(data))
    for i in range(68904):
        if (i % 72) < 30 or (i % 72) > 47:
            data_mask_row[i] = True
    data_mask = np.zeros(data.shape)


    data_mask[:, 0] = data_mask_row
    data_mask[:, 1] = data_mask_row
    data_mask[:, 2] = data_mask_row
    data_mask[:, 9] = data_mask_row
    data_mask[:, 10] = data_mask_row
    data_mask[:, 11] = data_mask_row


    dataframe = pp.DataFrame(data, mask=data_mask)
    datatime = np.arange(len(data))

    # tp.plot_timeseries(data, datatime, var_names, use_mask=True,
    #                    mask=data_mask, grey_masked_samples='data')

    parcorr = ParCorr(significance='analytic', use_mask=True, mask_type='y')
    pcmci = PCMCI(dataframe=dataframe, cond_ind_test=parcorr,
                  var_names=var_names, verbosity=1)

    # correlations = pcmci.get_lagged_dependencies(tau_max=20)
    # lag_func_matrix = tp.plot_lagfuncs(val_matrix=correlations,
    #                                    setup_args={'var_names':var_names, 
    #                                    'x_base':5, 'y_base':.5})
    pcmci.verbosity = 1
    results = pcmci.run_pcmci(tau_max=6, tau_min=1, pc_alpha=0.01)

    # print("p-values")
    # print (results['p_matrix'].round(3))
    # print("MCI partial correlations")
    # print (results['val_matrix'].round(2))

    q_matrix = pcmci.get_corrected_pvalues(p_matrix=results['p_matrix'],
                                           fdr_method='fdr_bh')
    pcmci._print_significant_links(p_matrix=results['p_matrix'],
                                   q_matrix=q_matrix,
                                   val_matrix=results['val_matrix'],
                                   alpha_level=0.01)

    link_matrix = pcmci._return_significant_parents(pq_matrix=q_matrix,
                                                    val_matrix=results['val_matrix'],
                                                    alpha_level=0.01)['link_matrix']

    tp.plot_time_series_graph(val_matrix=results['val_matrix'],
                              link_matrix=link_matrix, var_names=var_names,
                              link_colorbar_label='MCI',)
    return results, link_matrix


def caus_gpdc(data, var_names):
    import numpy as np
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    import sklearn

    import tigramite
    from tigramite import data_processing as pp
    from tigramite import plotting as tp
    from tigramite.pcmci import PCMCI
    from tigramite.independence_tests import ParCorr, GPDC, CMIknn, CMIsymb
    from tigramite.models import LinearMediation, Prediction

    data_mask_row = np.zeros(len(data))
    for i in range(68904):
        if (i % 72) < 30 or (i % 72) > 47:
            data_mask_row[i] = True
    data_mask = np.zeros(data.shape)

    data_mask[:, 0] = data_mask_row
    data_mask[:, 1] = data_mask_row
    data_mask[:, 2] = data_mask_row
    data_mask[:, 9] = data_mask_row
    data_mask[:, 10] = data_mask_row
    data_mask[:, 11] = data_mask_row

    dataframe = pp.DataFrame(data, mask=data_mask)
    datatime = np.arange(len(data))

    # tp.plot_timeseries(data, datatime, var_names, use_mask=True,
    #                    mask=data_mask, grey_masked_samples='data')

    gpdc = GPDC(significance='analytic', gp_params=None, use_mask=True, mask_type='y')
    gpdc.generate_and_save_nulldists(sample_sizes=range(495, 501),
                                     null_dist_filename='dc_nulldists.npz')
    gpdc.null_dist_filename ='dc_nulldists.npz'
    pcmci_gpdc = PCMCI(dataframe=dataframe, cond_ind_test=gpdc,
                       var_names=var_names, verbosity=1)

    # correlations = pcmci.get_lagged_dependencies(tau_max=20)
    # lag_func_matrix = tp.plot_lagfuncs(val_matrix=correlations,
    #                                    setup_args={'var_names':var_names,
    #                                    'x_base':5, 'y_base':.5})

    results = pcmci_gpdc.run_pcmci(tau_max=6, tau_min=1, pc_alpha=0.01)

    # print("p-values")
    # print (results['p_matrix'].round(3))
    # print("MCI partial correlations")
    # print (results['val_matrix'].round(2))

    q_matrix = pcmci_gpdc.get_corrected_pvalues(p_matrix=results['p_matrix'],
                                           fdr_method='fdr_bh')
    pcmci_gpdc._print_significant_links(p_matrix=results['p_matrix'],
                                   q_matrix=q_matrix,
                                   val_matrix=results['val_matrix'],
                                   alpha_level=0.01)

    link_matrix = pcmci_gpdc._return_significant_parents(pq_matrix=q_matrix,
                                                    val_matrix=results['val_matrix'],
                                                    alpha_level=0.01)['link_matrix']

    tp.plot_time_series_graph(val_matrix=results['val_matrix'],
                              link_matrix=link_matrix, var_names=var_names,
                              link_colorbar_label='MCI',)
    return results, link_matrix


def animatesu(data, colorlimit, precip='no',
              cbarleg='cbarlg', b1='batch_521', b2='batch_522'):
    """
    Plots animation of input grid of quantity at sigma and lat coords

    Parameters
    ----------
    data: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from netcdfread import ncread
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')

    summer = data[b2] - data[b1]

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    summer, lon3 = shiftgrid(180., summer, lon, start=False)
    summer, lon3 = addcyclic(summer, lon3)
    meshlon, meshlat = np.meshgrid(lon3, lat)
    ctrs = np.linspace(-colorlimit, colorlimit, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    if precip == 'yes':
        my_cmap = my_cmap[::-1]
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap

    def updatefig(i):
        fig.clear()

        ax1 = fig.add_subplot(1, 1, 1)
        m = Basemap(projection='moll', lon_0=0, resolution='c')
        m.drawcoastlines()
        m.drawcountries()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        plot = m.contourf(x, y, summer[i], ctrs,
                          cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        # parallels = m.drawparallels(np.arange(-90., 91., 15.))
        # meridians = m.drawmeridians(np.arange(-180., 181., 30.))
        # m.drawparallels(parallels, labels=[True, True, True, True])
        # m.drawmeridians(meridians, labels=[True, True, True, True])

        cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.015])
        b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                         orientation='horizontal', extend='max')
        b.set_label(label=cbarleg, size=20, fontsize=20, fontname='Arial')
        cl = plt.getp(cbar_ax, 'xmajorticklabels')
        ax1.set_title('Day ' + str(i), y=1.08, fontsize=20)
        plt.setp(cl, fontname='Arial', fontsize=20)

        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig, len(summer))
    # anim.save('test.mp4',bitrate=10000)
    return anim















