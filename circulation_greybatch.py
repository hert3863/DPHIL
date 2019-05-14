#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 14 15:08:36 2017

@author: bakerh
"""

import numpy as np


def compare_var(items):
    '''
    Compares the control and perturbed ensembles
    and outputs a list of exp IDs that have completed
    for both ensembles, coupled with the patch number
    '''
    import glob
    import fnmatch
    both = {}
    exps = ['batch_521', 'batch_522']
    for exp in exps:
        # import lists of successful files for each exp
        t521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/atmos/item15202_monthly_mean/*')
        i521 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/' + exp + '/atmos/' + items + '/*')

        # turn lists into just lists of exp IDs
        for i, item in enumerate(t521):
            t521[i] = t521[i][98+2*(len('item15202_monthly_mean')-22):102+2*(len('item15202_monthly_mean')-22)]
        for i, item in enumerate(i521):
            i521[i] = i521[i][98+2*(len(items)-22):102+2*(len(items)-22)]

        b = []
        # compare lists and add to dictionary if both exist
        for i, item in enumerate(i521):
            if fnmatch.filter(t521, i521[i]) != []:
                b.append(i521[i])
        both[exp] = b
    return both


def importmeans(item, both, daily='no', wind='no'):
    from netCDF4 import Dataset
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

    outputv = {}
    outputf = {}
    exps = ['batch_521', 'batch_522']

    for x, exp in enumerate(exps):
        a = []
        for i in range(np.ma.size(both[exp])):
            a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                     '/atmos/item15202_monthly_mean/item15202_monthly_mean_' +
                     both[exp][i] + '_2090-01_2100-12.nc')
        vbar = np.zeros((2, np.ma.size(a), 144, 192))
        sdata = np.zeros((3*11*np.ma.size(a), 144, 192))
        wdata = np.zeros((3*11*np.ma.size(a), 144, 192))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata[33*e:33*(e+1)] = nc_fid.variables['item15202_monthly_mean'
                                                    ][sindices, 2, :]
            wdata[33*e:33*(e+1)] = nc_fid.variables['item15202_monthly_mean'
                                                    ][windices, 2, :]
            print('Done ' + str('item15202_monthly_mean') + ' ' +
                  str(exp) + ' ' + str(e+1))
        for y in range(np.ma.size(a)):
            vbar[0, y, :] = np.mean(wdata[33*y:33*(y+1), :], axis=0)
            vbar[1, y, :] = np.mean(sdata[33*y:33*(y+1), :], axis=0)
        outputv[exp] = vbar

    if wind == 'yes':
        le = 144
        lev = 2
    else:
        le = 145
        lev = 0
    if daily == 'yes':
        for x, exp in enumerate(exps):
            a = []
            for i in range(np.ma.size(both[exp])):
                a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                         '/atmos/' + item + '/' + item + '_' +
                         both[exp][i] + '_2090-01_2100-12.nc')
            rbar = np.zeros((2, np.ma.size(a), le, 192))
            sdatad = np.zeros((90*11*np.ma.size(a), le, 192))
            wdatad = np.zeros((90*11*np.ma.size(a), le, 192))

            for e, d in enumerate(a):
                nc_fid = Dataset(d, 'r')
                sdatad[990*e:990*(e+1)] = nc_fid.variables[item
                                                           ][sindicesd, lev, :]
                wdatad[990*e:990*(e+1)] = nc_fid.variables[item
                                                           ][windicesd, lev, :]
                print('Done ' + str(item) + ' ' +
                      str(exp) + ' ' + str(e+1))
            for y in range(np.ma.size(a)):
                rbar[0, y, :] = np.mean(wdatad[990*y:990*(y+1), :], axis=0)
                rbar[1, y, :] = np.mean(sdatad[990*y:990*(y+1), :], axis=0)

            outputf[exp] = rbar * 86400
    else:
        for x, exp in enumerate(exps):
            a = []
            for i in range(np.ma.size(both[exp])):
                a.append('/network/aopp/hera/mad/bakerh/HAPPI/' + exp +
                         '/atmos/' + item + '/' + item + '_' +
                         both[exp][i] + '_2090-01_2100-12.nc')
            fbar = np.zeros((2, np.ma.size(a), le, 192))
            sdata = np.zeros((3*11*np.ma.size(a), le, 192))
            wdata = np.zeros((3*11*np.ma.size(a), le, 192))
            for e, d in enumerate(a):
                nc_fid = Dataset(d, 'r')
                sdata[33*e:33*(e+1)] = nc_fid.variables[item
                                                        ][sindices, lev, :]
                wdata[33*e:33*(e+1)] = nc_fid.variables[item
                                                        ][windices, lev, :]
                print('Done ' + str(item) + ' ' +
                      str(exp) + ' ' + str(e+1))
            for y in range(np.ma.size(a)):
                fbar[0, y, :] = np.mean(wdata[33*y:33*(y+1), :], axis=0)
                fbar[1, y, :] = np.mean(sdata[33*y:33*(y+1), :], axis=0)
            outputf[exp] = fbar
    return outputv, outputf


def pcd(z500, r):
    from scipy import stats
    z500_lower = np.mean(z500['batch_521'][1], axis=0)
    z500_anom = z500['batch_522'][1] - z500_lower
    # z500_anom = shiftlat(z500_anom)
    z500_anom = (z500_anom - np.mean(z500_anom, axis=0))/np.std(z500_anom, axis=0)

    r_lower = np.mean(r['batch_521'][1], axis=0)
    r_anom = r['batch_522'][1] - r_lower
    r_pcd = np.mean(r_anom[:, 56:65, 96:134],
                    axis=(1, 2)) - np.mean(r_anom[:, 52:65, 147:158],
                                           axis=(1, 2))
    beta = np.zeros((145, 192))
    for i in range(145):
        for j in range(192):
            beta[i, j] = stats.linregress(r_pcd[:], z500_anom[:, i, j])[0]
    return beta


def shiftlat(grid):
    latold = np.arange(89.375, -90.625, -1.25)
    latnew = np.arange(90, -91.25, -1.25)
    regrid = np.zeros((np.ma.size(grid, axis=0), 145, 192))
    for i in range(143):
        regrid[:, i+1, :] = ((grid[:, i, :]*np.cos(latold[i]*np.pi/180) +
                              grid[:, i+1, :]*np.cos(latold[i+1]*np.pi/180)) /
                             (2*np.cos(latnew[i+1]*np.pi/180)))
    return regrid


def v_reg(v, r):
    from scipy import stats
    v_lower = np.mean(v['batch_521'][1], axis=0)
    v_anom = v['batch_522'][1] - v_lower
    v_anom = shiftlat(v_anom)
    r_lower = np.mean(r['batch_521'][1], axis=0)
    r_anom = r['batch_522'][1] - r_lower
    r_anom = shiftlat(r_anom)
    # r_anom = (r_anom - np.mean(r_anom, axis=0))/np.std(r_anom, axis=0)
    v_ind = v_anom[:, 30, 0]
    v_ind = (v_ind - np.mean(v_ind, axis=0))/np.std(v_ind, axis=0)
    beta = np.zeros((145, 192))
    for i in range(145):
        for j in range(192):
            beta[i, j] = stats.linregress(v_ind[:], r_anom[:, i, j])[0]
    return beta








































