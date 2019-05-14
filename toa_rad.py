# -*- coding: utf-8 -*-
"""
Created on Tue Oct  3 14:43:14 2017

@author: bakerh
"""


import numpy as np
from netcdfread import *
import glob


def radbal():
    a = [sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/AMIP/rsdt/*'), key=str.lower),
         sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/AMIP/rsut/*'), key=str.lower),
         sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/AMIP/rlut/*'), key=str.lower)]
    N = np.zeros((np.ma.size(a[0])))
    name = []
    for i, item in enumerate(a[0]):
        print(a[0][i][56:-29])
        lon = ncread(item, 'lon')
        lat = ncread(item, 'lat')
        s = 0
        e = 360
        if a[0][i][56:-29] == 'ACCESS1-3':
            s = 12
            e = 372
        if a[0][i][56:-29] == 'CanAM4':
            s = -372
            e = -12
        if a[0][i][56:-29] == 'CCSM4':
            s = 0
            e = -24
        if a[0][i][56:-29] == 'CSIRO-Mk3-6-0':
            s = 0
            e = -12
        if a[0][i][56:-29] == 'FGOALS-g2':
            s = 0
            e = -12
        if a[0][i][56:-29] == 'GISS-E2-R':
            s = -384
            e = -24
        if a[0][i][56:-29] == 'HadGEM2-A':
            s = -360
            e = 364
        if a[0][i][56:-29] == 'IPSL-CM5A-LR':
            s = 0
            e = -12
        if a[0][i][56:-29] == 'IPSL-CM5A-MR':
            s = -372
            e = -12
        rsdt = ncread(a[0][i], 'rsdt')[s:e]
        rsut = ncread(a[1][i], 'rsut')[s:e]
        rlut = ncread(a[2][i], 'rlut')[s:e]
        meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
        meshlat[:, :] = lat
        meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
        rsdt *= meshlatweight
        rsdt = np.mean(rsdt) / np.mean(meshlatweight)
        rsut *= meshlatweight
        rsut = np.mean(rsut) / np.mean(meshlatweight)
        rlut *= meshlatweight
        rlut = np.mean(rlut) / np.mean(meshlatweight)
        N[i] = rsdt-rsut-rlut
        name.append(a[0][i][56:-29])
    return N, name


def scttr(n_amip, n_had, n_cam, names):
    from matplotlib import pyplot as plt
    from adjustText import adjust_text
    hfont = {'fontname': 'Arial'}
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = [1, 2, 3]
    labels = ['AMIP', 'HadAM3P', 'CAM4']
    names_had = ['Present', 'B.E.', 'Low', 'High']
    names_cam = ['Present', 'B.E.', 'Low']
    plt.scatter(np.ones((24)), n_amip, color='red',s=50)
    plt.scatter(1, n_amip.mean(), color='black',s=50)
    plt.scatter([2, 2, 2, 2], n_had, color='blue',s=50)
    plt.scatter([3, 3, 3], n_cam, color='green',s=50)
    #color = ['b', 'k', 'r']
    #pp = []
    ax1.set_xlim([0, 4])
    ax1.set_ylim([-2.5, 6])
    ax1.set_ylabel('Net TOA downward radiation (Wm$\mathregular{^{-2}}$)',
                   fontsize=20, **hfont)
    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation='vertical', **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(4, linewidth=2.5, color='k')
    ax1.axhline(0, linewidth=1.5, color='k')
    ax1.axhline(-2.5, linewidth=2.5, color='k')
    ax1.axhline(6, linewidth=2.5, color='k')
    '''
    for label, x1, y1 in zip(names, np.ones((24)), n_amip):
        ax1.annotate(label, xy=(x1, y1), xytext=(-.20, -.1),
                     textcoords='offset points', ha='right', va='bottom',
                     fontsize=12)
    '''
    ax1.annotate('Present', (1.6, n_had[0]+.05))
    ax1.annotate('B.E.', (1.7, n_had[1]-.1))
    ax1.annotate('Low', (1.7, n_had[2]-.05))
    ax1.annotate('High', (1.7, n_had[3]-.05))
    ax1.annotate('Present', (2.6, n_cam[0]-.05))
    ax1.annotate('B.E.', (2.7, n_cam[1]-.05))
    ax1.annotate('Low', (2.7, n_cam[2]-.05))
    '''
    for label, x1, y1 in zip(names_had, [2, 2, 2, 2], n_had):
        ax1.annotate(label, xy=(x1, y1), xytext=(-.20, -.1),
                     #textcoords='offset points', ha='right', va='bottom',
                      fontsize=12)
    for label, x1, y1 in zip(names_cam, [3, 3, 3], n_cam):
        ax1.annotate(label, xy=(x1, y1), xytext=(-.5, -.5),
                     #textcoords='offset points', ha='right', va='bottom',
                     fontsize=12)
    #ax1.legend(handles=pp, labels=['Low CO$\mathregular{_{2}}$',
       #                            'Mean CO$\mathregular{_{2}}$',
        #                           'High CO$\mathregular{_{2}}$'],
       #        loc=3, scatterpoints=1, prop={'family': 'Arial'},
         #      bbox_to_anchor=(.885, .86, 1., .102))
    '''
    plt.subplots_adjust(hspace=.3, wspace=0.05, top=.97, bottom=0.2, left=.3,
                        right=.9)


















