# -*- coding: utf-8 -*-
"""
Created on Fri Sep  8 17:11:18 2017

@author: bakerh
"""

import numpy as np


def fig4_multi1(tmeans, conversion, rcpplume, alpha, beta, cvs,
               tmeans1, alpha1, beta1, cvs1,
               tmeans2, alpha2, beta2, cvs2,
               title=''):
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
    # co2 = np.array([390.4, 423.4, 395.8, 550.0])
    tmean = np.zeros((4))
    tmean[0] = np.mean(tmeans['batch_518'])
    tmean[1] = np.mean(tmeans['batch_520'])
    tmean[2] = np.mean(tmeans['batch_521'])
    tmean[3] = np.mean(tmeans['batch_522'])

    tmean1 = np.zeros((4))
    tmean1[0] = np.mean(tmeans1['batch_518'])
    tmean1[1] = np.mean(tmeans1['batch_520'])+0.05
    tmean1[2] = np.mean(tmeans1['batch_521'])
    tmean1[3] = np.mean(tmeans1['batch_522'])

    tmean2 = np.zeros((4))
    tmean2[0] = np.mean(tmeans2['batch_518'])
    tmean2[1] = np.mean(tmeans2['batch_520'])-.05
    tmean2[2] = np.mean(tmeans2['batch_521'])
    tmean2[3] = np.mean(tmeans2['batch_522'])

    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv['d_f'][1] +
                       (x*alpha/beta**2)**2*cv['d_t'][0] -
                       2*(x**2*alpha*cv['d_t'][1]**2/beta**3))
        return sigT

    # co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301

    ax1.set_ylim([0, 1.5])
    hfont = {'fontname': 'Arial'}
    ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2015 \
($^{\circ}$C)', fontsize=20, **hfont)
    ax1.set_xlabel('Cumulative total anthropogenic CO$\mathregular{_{2}}$ \
emissions from 2011 (Gt C)', fontsize=20, **hfont)
    ax1.set_xticks(np.arange(0, 1100, 100))
    ax1.set_xlim([0, 650])
    # ax1.plot(conversion[:, 2], 0.01951515*(conversion[:, 0]-2010))
    # ax1.plot(conversion[156:162, 2], test1[4:])
    x11 = rcpplume[0, 20]
    y11 = rcpplume[1, 20]
    x22 = rcpplume[0, 3]
    y22 = rcpplume[1, 3]
    g1 = (rcpplume[1, 19]-rcpplume[1, 20])/(rcpplume[0, 19]-rcpplume[0, 20])
    g2 = (rcpplume[1, 4]-rcpplume[1, 3])/(rcpplume[0, 4]-rcpplume[0, 3])
    ax1.axhline(tmean[1]-tmean[0], linestyle='--', color='#0708EC', linewidth=2,
                dashes=(10, 10))
    ax1.fill_between(rcpplume[0, :], rcpplume[1, :], facecolor='#F5D9CE',
                     edgecolor='#eebdaa', alpha=1)
    x = np.arange(0, 1010, .01)

    # first ensemble
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
    spr = ax1.axvspan(ubl, ubu, color='#0708EC', alpha=0.25)
    ax1.axvline(oub, linestyle='--', color='#0708EC', linewidth=2)
    ax1.axvline(lb, linestyle='--', color='#0708EC', linewidth=2)
    ubp = ax1.axvline(ub, linestyle='--', color='#0708EC', linewidth=2)
    mns = ax1.scatter(rcp_conc, tmean[1]-tmean[0], color='#5EB3FA',
                      marker='o', s=100, lw=3, label='RCP2.6 2090s')
    # second ensemble
    ax1.axhline(tmean1[1]-tmean1[0], linestyle='--', color='#ff6600', linewidth=2,
                dashes=(10, 10))
    y_1 = (-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0])
    cns1, = ax1.plot(x, y_1, linewidth=2, color='#ff6600')
    y1_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) -
            1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
    y2_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) +
            1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
    ax1.fill_between(x, y1_1, y2_1, facecolor='#ff6600', edgecolor='#ff6600',
                     alpha=.25)
    oub1 = (tmean1[1]-tmean1[0] - y22) / g2 + x22
    lb1 = (tmean1[1]-tmean1[0] - y11) / g1 + x11
    ub1 = x[np.argwhere(np.diff(np.sign(y_1 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubu1 = x[np.argwhere(np.diff(np.sign(y2_1 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubl1 = x[np.argwhere(np.diff(np.sign(y1_1 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    spr1 = ax1.axvspan(ubl1, ubu1, color='#ff6600', alpha=0.25)
    ax1.axvline(oub1, linestyle='--', color='#ff6600', linewidth=2)
    ax1.axvline(lb1, linestyle='--', color='#ff6600', linewidth=2)
    ubp1 = ax1.axvline(ub1, linestyle='--', color='#ff6600', linewidth=2)
    mns1 = ax1.scatter(rcp_conc, tmean1[1]-tmean1[0], color='#5EB3FA',
                       marker='o', s=100, lw=3, label='RCP2.6 2090s')
    # third ensemble

    ax1.axhline(tmean2[1]-tmean2[0], linestyle='--', color='#FF0400', linewidth=2,
                dashes=(10, 10))
    y_2 = (-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0])
    cns2, = ax1.plot(x, y_2, linewidth=2, color='#FF0400')
    y1_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) -
            1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
    y2_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) +
            1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
    ax1.fill_between(x, y1_2, y2_2, facecolor='#FF0400', edgecolor='#FF0400',
                     alpha=.25)
    oub2 = (tmean2[1]-tmean2[0] - y22) / g2 + x22
    lb2 = (tmean2[1]-tmean2[0] - y11) / g1 + x11
    ub2 = x[np.argwhere(np.diff(np.sign(y_2 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubu2 = x[np.argwhere(np.diff(np.sign(y2_2 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubl2 = x[np.argwhere(np.diff(np.sign(y1_2 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    spr2 = ax1.axvspan(ubl2, ubu2, color='#FF0400', alpha=0.25)
    ax1.axvline(oub2, linestyle='--', color='#FF0400', linewidth=2)
    ax1.axvline(lb2, linestyle='--', color='#FF0400', linewidth=2)
    ubp2 = ax1.axvline(ub2, linestyle='--', color='#FF0400', linewidth=2)
    mns2 = ax1.scatter(rcp_conc, tmean2[1]-tmean2[0], color='#5EB3FA',
                       marker='o', s=100, lw=3, label='RCP2.6 2090s')

    # as normal
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
                       'Constant TX90p',
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


def fig4_multi2(tmean, conversion, rcpplume, alpha, beta, cvs,
                tmean1, alpha1, beta1, cvs1,
                tmean2, alpha2, beta2, cvs2,
                title=''):
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
    # co2 = np.array([390.4, 423.4, 395.8, 550.0])
    '''
    tmean = np.zeros((4))
    tmean[0] = np.mean(tmeans['batch_518'])
    tmean[1] = np.mean(tmeans['batch_520'])
    tmean[2] = np.mean(tmeans['batch_521'])
    tmean[3] = np.mean(tmeans['batch_522'])

    tmean1 = np.zeros((3))
    tmean1[0] = np.mean(tmeans1['present'])
    tmean1[1] = np.mean(tmeans1['plus15'])
    tmean1[2] = np.mean(tmeans1['plus15_lower'])
    tmean1[3] = np.mean(tmeans1['plus15_higher'])

    tmean2 = np.zeros((4))
    tmean2[0] = np.mean(tmeans2['batch_518'])
    tmean2[1] = np.mean(tmeans2['batch_520'])+.05
    tmean2[2] = np.mean(tmeans2['batch_521'])
    tmean2[3] = np.mean(tmeans2['batch_522'])
    '''
    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv[1, 1] +
                       (x*alpha/beta**2)**2*cv[0, 0] -
                       2*(x**2*alpha*cv[1, 0]**2/beta**3))
        return sigT

    # co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301

    ax1.set_ylim([0, 1.4])
    hfont = {'fontname': 'Arial'}
    ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2015 \
($^{\circ}$C)', fontsize=20, **hfont)
    ax1.set_xlabel('Cumulative total anthropogenic CO$\mathregular{_{2}}$ \
emissions from 2011 (Gt C)', fontsize=20, **hfont)
    ax1.set_xticks(np.arange(0, 1100, 100))
    ax1.set_xlim([0, 650])
    # ax1.plot(conversion[:, 2], 0.01951515*(conversion[:, 0]-2010))
    # ax1.plot(conversion[156:162, 2], test1[4:])
    x11 = rcpplume[0, 20]
    y11 = rcpplume[1, 20]
    x22 = rcpplume[0, 3]
    y22 = rcpplume[1, 3]
    g1 = (rcpplume[1, 19]-rcpplume[1, 20])/(rcpplume[0, 19]-rcpplume[0, 20])
    g2 = (rcpplume[1, 4]-rcpplume[1, 3])/(rcpplume[0, 4]-rcpplume[0, 3])
    ax1.axhline((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3,
                linestyle='--', color='#091034',
                linewidth=2, dashes=(10, 10))
    # bls = ax1.axhspan(tmean1[1]-tmean1[0], tmean2[1]-tmean2[0],
    #                   color='#091034', alpha=0.25)
    ax1.fill_between(rcpplume[0, :], rcpplume[1, :], facecolor='#F7E3D4',
                     edgecolor='#F7E3D4', alpha=1)
    x = np.arange(0, 1010, .01)

    # first ensemble
    y = (-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0])
    cns, = ax1.plot(x, y, linewidth=2, color='#193DF0')
    y1 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) -
          1.6445*sigt(gtc2df(x), cvs, alpha, beta))
    y2 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) +
          1.6445*sigt(gtc2df(x), cvs, alpha, beta))
    ax1.fill_between(x, y1, y2, facecolor='#193DF0', edgecolor='#193DF0',
                     alpha=.25)
    oub = ((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3 -
           y22) / g2 + x22
    lb = ((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3 -
          y11) / g1 + x11
    ub = x[np.argwhere(np.diff(np.sign(y -
                               g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubu = x[np.argwhere(np.diff(np.sign(y2 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubl = x[np.argwhere(np.diff(np.sign(y1 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]

    ax1.axvline(oub, linestyle='--', color='#FD3D12', linewidth=2)
    ax1.axvline(lb, linestyle='--', color='#FD3D12', linewidth=2)

    mns = ax1.scatter(rcp_conc, tmean[1]-tmean[0], color='#193DF0',
                      marker='o', s=100, lw=3, label='RCP2.6 2090s')
    # second ensemble
    y_1 = (-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0])
    cns1, = ax1.plot(x, y_1, linewidth=2, color='#2369B8')
    y1_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) -
            1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
    y2_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) +
            1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
    ax1.fill_between(x, y1_1, y2_1, facecolor='#2369B8', edgecolor='#2369B8',
                     alpha=.25)
    # oub1 = (tmean1[1]-tmean1[0] - y22) / g2 + x22
    # lb1 = (tmean1[1]-tmean1[0] - y11) / g1 + x11
    ub1 = x[np.argwhere(np.diff(np.sign(y_1 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubu1 = x[np.argwhere(np.diff(np.sign(y2_1 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubl1 = x[np.argwhere(np.diff(np.sign(y1_1 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]

    mns1 = ax1.scatter(rcp_conc, tmean1[1]-tmean1[0], color='#2369B8',
                       marker='o', s=100, lw=3, label='RCP2.6 2090s')
    # third ensemble
    y_2 = (-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0])
    cns2, = ax1.plot(x, y_2, linewidth=2, color='#5918C9')
    y1_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) -
            1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
    y2_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) +
            1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
    ax1.fill_between(x, y1_2, y2_2, facecolor='#5918C9', edgecolor='#5918C9',
                     alpha=.25)
    # oub2 = (tmean2[1]-tmean2[0] - y22) / g2 + x22
    # lb2 = (tmean2[1]-tmean2[0] - y11) / g1 + x11
    ub2 = x[np.argwhere(np.diff(np.sign(y_2 -
                                g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubu2 = x[np.argwhere(np.diff(np.sign(y2_2 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    ubl2 = x[np.argwhere(np.diff(np.sign(y1_2 -
                                 g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    spr = ax1.axvspan(np.min([ubl, ubl1, ubl2]), np.max([ubu, ubu1, ubu2]),
                      color='#FBAA09', alpha=0.25)
    #spr = ax1.axvspan((ub+ub1+ub2)/3-np.sqrt([(ubl-ub)**2+(ubl1-ub1)**2+(ubl2-ub2)**2])/3, (ub+ub1+ub2)/3+np.sqrt([(ubu-ub)**2+(ubu1-ub1)**2+(ubu2-ub2)**2])/3,
     #                 color='#FBAA09', alpha=0.25)
    #spr = ax1.axvspan((ub+ub1+ub2)/3-1.6445*np.std([ub,ub1,ub2]), (ub+ub1+ub2)/3+1.6445*np.std([ub,ub1,ub2]),
     #                 color='#FBAA09', alpha=0.25)
    # rls = ax1.axvspan(oub2, oub1, color='#FD3D12', alpha=0.25)
    # ax1.axvspan(lb2, lb1, color='#FD3D12', alpha=0.25)
    ubp = ax1.axvline((ub+ub1+ub2)/3, linestyle='--', color='#FBAA09',
                      linewidth=2)
    mns2 = ax1.scatter(rcp_conc, tmean2[1]-tmean2[0], color='#5918C9',
                       marker='o', s=100, lw=3, label='RCP2.6 2090s')

    # as normal
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(ax1.get_xticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    # ax2.xaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(1000, linewidth=2.5, color='k')
    ax1.axhline(0, linewidth=2.5, color='k')
    ax1.axhline(1.4, linewidth=2.5, color='k')

    pink_patch = mpatches.Patch(color='#F7E3D4', label='RCP range')
    red_line = mlines.Line2D([], [], color='#FD3D12', linestyle='--',
                             linewidth=2, label='Carbon budget bound')
    black_line = mlines.Line2D([], [], color='#091034', linestyle='--',
                               linewidth=2,
                               dashes=(10, 10), label='1.5$^{\circ}$C warming')
    b_spr = mpatches.Patch(color='#193DF0', alpha=0.25, linewidth=0)
    b_spr1 = mpatches.Patch(color='#2369B8', alpha=0.25, linewidth=0)
    b_spr2 = mpatches.Patch(color='#5918C9', alpha=0.25, linewidth=0)
    ax1.legend(handles=[pink_patch, black_line, mns, (cns, b_spr), mns1,
                        (cns1, b_spr1), mns2, (cns2, b_spr2),
                        red_line, (ubp, spr)],
               labels=['RCP range', '1.5$^{\circ}$C warming',
                       'RCP2.6 2090s HadAM3P', 'Constant TX90p HadAM3P',
                       'RCP2.6 2090s MIROC5', 'Constant TX90p MIROC5',
                       'RCP2.6 2090s CAM4', 'Constant TX90p CAM4',
                       'Carbon budget bound', 'New upper bound'],
               handlelength=3,
               loc=3, scatterpoints=1, prop={'family': 'Arial'},
               bbox_to_anchor=(.675, .65, 2, 1), frameon=True)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.985)
    plt.title(title, y=1.01, fontsize=20, **hfont)
    # plt.subplots_adjust(hspace=0, wspace=0.05, top=.95, bottom=0.1, left=.05,
    #                   right=.97)
    un = oub - lb
    ub_m = (ub+ub1+ub2)/3
    ubl_l = np.min([ubl, ubl1, ubl2])
    ubu_u = np.max([ubu, ubu1, ubu2])
    # ubl_l = (ub+ub1+ub2)/3-np.sqrt([(ubl-ub)**2+(ubl1-ub1)**2+(ubl2-ub2)**2])/3
    # ubu_u = (ub+ub1+ub2)/3+np.sqrt([(ubu-ub)**2+(ubu1-ub1)**2+(ubu2-ub2)**2])/3
    p_ch = np.array([float(ub_m)-oub, float(ub_m)-float(ubl_l),
                     float(ubu_u)-float(ub_m)])*100/un

    return oub, float(ub_m), float(ubl_l), float(ubu_u), p_ch


def fig4_multi3(tmean, conversion, rcpplume, a, b, c,
                tmean1, a1, b1, c1,
                tmean2, a2, b2, c2,
                title=''):
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
    fs = 12
    fig, ax = plt.subplots(2, 2, figsize=(13, 13),
                           facecolor='w', edgecolor='k', linewidth=2)
    fig.delaxes(ax[1, 1])
    # co2 = np.array([390.4, 423.4, 395.8, 550.0])
    '''
    tmean = np.zeros((4))
    tmean[0] = np.mean(tmeans['batch_518'])
    tmean[1] = np.mean(tmeans['batch_520'])
    tmean[2] = np.mean(tmeans['batch_521'])
    tmean[3] = np.mean(tmeans['batch_522'])

    tmean1 = np.zeros((3))
    tmean1[0] = np.mean(tmeans1['present'])
    tmean1[1] = np.mean(tmeans1['plus15'])
    tmean1[2] = np.mean(tmeans1['plus15_lower'])
    tmean1[3] = np.mean(tmeans1['plus15_higher'])

    tmean2 = np.zeros((4))
    tmean2[0] = np.mean(tmeans2['batch_518'])
    tmean2[1] = np.mean(tmeans2['batch_520'])+.05
    tmean2[2] = np.mean(tmeans2['batch_521'])
    tmean2[3] = np.mean(tmeans2['batch_522'])
    '''
    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv[1, 1] +
                       (x*alpha/beta**2)**2*cv[0, 0] -
                       2*(x**2*alpha*cv[1, 0]**2/beta**3))
        return sigT

    # co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301
    var = ['R95p', 'TX90p', 'WBGT95p']
    let = ['c', 'a', 'b']
    for i in range(3):
        ax1 = ax[int(np.floor(i/2)), int(np.floor(1/(i+1)))]
        alpha = a[2-i]
        beta = b[2-i]
        cvs = c[2-i]
        alpha1 = a1[2-i]
        beta1 = b1[2-i]
        cvs1 = c1[2-i]
        alpha2 = a2[2-i]
        beta2 = b2[2-i]
        cvs2 = c2[2-i]
        v = var[2-i]
        l = let[2-i]
        ax1.set_ylim([0, 1.4])
        hfont = {'fontname': 'Arial'}
        ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2015 \
($^{\circ}$C)', fontsize=fs, **hfont)
        ax1.set_xlabel('Cumulative total anthropogenic CO$\mathregular{_{2}}$ \
emissions from 2011 (Gt C)', fontsize=fs, **hfont)
        ax1.set_xticks(np.arange(0, 1100, 100))
        ax1.set_xlim([0, 650])
        # ax1.plot(conversion[:, 2], 0.01951515*(conversion[:, 0]-2010))
        # ax1.plot(conversion[156:162, 2], test1[4:])
        x11 = rcpplume[0, 20]
        y11 = rcpplume[1, 20]
        x22 = rcpplume[0, 3]
        y22 = rcpplume[1, 3]
        g1 = (rcpplume[1, 19]-rcpplume[1, 20])/(rcpplume[0, 19]-rcpplume[0, 20])
        g2 = (rcpplume[1, 4]-rcpplume[1, 3])/(rcpplume[0, 4]-rcpplume[0, 3])
        ax1.axhline((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3,
                    linestyle='--', color='#091034',
                    linewidth=2)
        # bls = ax1.axhspan(tmean1[1]-tmean1[0], tmean2[1]-tmean2[0],
        #                   color='#091034', alpha=0.25)
        ax1.fill_between(rcpplume[0, :], rcpplume[1, :], facecolor='#F7E3D4',
                         edgecolor='#F7E3D4', alpha=1)
        x = np.arange(0, 1010, .01)
    
        # first ensemble
        y = (-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0])
        cns, = ax1.plot(x, y, linewidth=2, color='#193DF0')
        y1 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) -
              1.6445*sigt(gtc2df(x), cvs, alpha, beta))
        y2 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) +
              1.6445*sigt(gtc2df(x), cvs, alpha, beta))
        ax1.fill_between(x, y1, y2, facecolor='#193DF0', edgecolor='#193DF0',
                         alpha=.25)
        oub = ((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3 -
               y22) / g2 + x22
        lb = ((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3 -
              y11) / g1 + x11
        ub = x[np.argwhere(np.diff(np.sign(y -
                                   g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu = x[np.argwhere(np.diff(np.sign(y2 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl = x[np.argwhere(np.diff(np.sign(y1 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    
        ax1.axvline(oub, linestyle='--', color='#FD3D12', linewidth=2)
        ax1.axvline(lb, linestyle='--', color='#FD3D12', linewidth=2)
    
        mns = ax1.scatter(rcp_conc, tmean[1]-tmean[0], color='#193DF0',
                          marker='o', s=100, lw=3, label='RCP2.6 2090s')
        # second ensemble
        y_1 = (-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0])
        cns1, = ax1.plot(x, y_1, linewidth=2, color='#2369B8')
        y1_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) -
                1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
        y2_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) +
                1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
        ax1.fill_between(x, y1_1, y2_1, facecolor='#2369B8', edgecolor='#2369B8',
                         alpha=.25)
        # oub1 = (tmean1[1]-tmean1[0] - y22) / g2 + x22
        # lb1 = (tmean1[1]-tmean1[0] - y11) / g1 + x11
        ub1 = x[np.argwhere(np.diff(np.sign(y_1 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu1 = x[np.argwhere(np.diff(np.sign(y2_1 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl1 = x[np.argwhere(np.diff(np.sign(y1_1 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    
        mns1 = ax1.scatter(rcp_conc, tmean1[1]-tmean1[0], color='#2369B8',
                           marker='o', s=100, lw=3, label='RCP2.6 2090s')
        # third ensemble
        y_2 = (-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0])
        cns2, = ax1.plot(x, y_2, linewidth=2, color='#5918C9')
        y1_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) -
                1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
        y2_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) +
                1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
        ax1.fill_between(x, y1_2, y2_2, facecolor='#5918C9', edgecolor='#5918C9',
                         alpha=.25)
        # oub2 = (tmean2[1]-tmean2[0] - y22) / g2 + x22
        # lb2 = (tmean2[1]-tmean2[0] - y11) / g1 + x11
        ub2 = x[np.argwhere(np.diff(np.sign(y_2 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu2 = x[np.argwhere(np.diff(np.sign(y2_2 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl2 = x[np.argwhere(np.diff(np.sign(y1_2 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        spr = ax1.axvspan(np.min([ubl, ubl1, ubl2]), np.max([ubu, ubu1, ubu2]),
                          color='#FBAA09', alpha=0.25)
        #spr = ax1.axvspan((ub+ub1+ub2)/3-np.sqrt([(ubl-ub)**2+(ubl1-ub1)**2+(ubl2-ub2)**2])/3, (ub+ub1+ub2)/3+np.sqrt([(ubu-ub)**2+(ubu1-ub1)**2+(ubu2-ub2)**2])/3,
         #                 color='#FBAA09', alpha=0.25)
        #spr = ax1.axvspan((ub+ub1+ub2)/3-1.6445*np.std([ub,ub1,ub2]), (ub+ub1+ub2)/3+1.6445*np.std([ub,ub1,ub2]),
         #                 color='#FBAA09', alpha=0.25)
        # rls = ax1.axvspan(oub2, oub1, color='#FD3D12', alpha=0.25)
        # ax1.axvspan(lb2, lb1, color='#FD3D12', alpha=0.25)
        ubp = ax1.axvline((ub+ub1+ub2)/3, linestyle='--', color='#FBAA09',
                          linewidth=2)
        mns2 = ax1.scatter(rcp_conc, tmean2[1]-tmean2[0], color='#5918C9',
                           marker='o', s=100, lw=3, label='RCP2.6 2090s')
    
        # as normal
        ax1.set_yticklabels(ax1.get_yticks(), fontsize=fs, **hfont)
        ax1.set_xticklabels(ax1.get_xticks(), fontsize=fs, **hfont)
        ax1.xaxis.set_tick_params(width=2, length=7.5)
        ax1.yaxis.set_tick_params(width=2, length=7.5)
        # ax2.xaxis.set_tick_params(width=2, length=7.5)
        ax1.axvline(0, linewidth=2.5, color='k')
        ax1.axvline(1000, linewidth=2.5, color='k')
        ax1.axhline(0, linewidth=2.5, color='k')
        ax1.axhline(1.4, linewidth=2.5, color='k')
    
        pink_patch = mpatches.Patch(color='#F7E3D4', label='RCP range')
        red_line = mlines.Line2D([], [], color='#FD3D12', linestyle='--',
                                 linewidth=2, label='Carbon budget bound')
        black_line = mlines.Line2D([], [], color='#091034', linestyle='--',
                                   linewidth=2, label='1.5$^{\circ}$C warming')
        b_spr = mpatches.Patch(color='#193DF0', alpha=0.25, linewidth=0)
        b_spr1 = mpatches.Patch(color='#2369B8', alpha=0.25, linewidth=0)
        b_spr2 = mpatches.Patch(color='#5918C9', alpha=0.25, linewidth=0)
        ax1.legend(handles=[pink_patch, black_line, mns, (cns, b_spr), mns1,
                            (cns1, b_spr1), mns2, (cns2, b_spr2),
                            red_line, (ubp, spr)],
                   labels=['RCP range', '1.5$^{\circ}$C warming',
                           'RCP2.6 2090s HadAM3P', 'Constant '+v+' HadAM3P',
                           'RCP2.6 2090s MIROC5', 'Constant '+v+' MIROC5',
                           'RCP2.6 2090s CAM4', 'Constant '+v+' CAM4',
                           'Carbon budget bound', 'New upper bound'],
                   handlelength=3,
                   loc=1, scatterpoints=1, prop={'family': 'Arial', 'size':fs-4},
                   bbox_to_anchor=(.995, .995), frameon=True, framealpha=1)
        ax1.set_title(l, loc='left', fontname='Arial', fontsize=fs+2,
                      fontweight='bold', x=-.15, y=1)
    plt.savefig('/home/bakerh/Downloads/fig3.png', dpi=300)
    #plt.title(title, y=1.01, fontsize=20, **hfont)
    # plt.subplots_adjust(hspace=0, wspace=0.05, top=.95, bottom=0.1, left=.05,
    #                   right=.97)
    #un = oub - lb
    #ub_m = (ub+ub1+ub2)/3
    #ubl_l = np.min([ubl, ubl1, ubl2])
    #ubu_u = np.max([ubu, ubu1, ubu2])
    ## ubl_l = (ub+ub1+ub2)/3-np.sqrt([(ubl-ub)**2+(ubl1-ub1)**2+(ubl2-ub2)**2])/3
    ## ubu_u = (ub+ub1+ub2)/3+np.sqrt([(ubu-ub)**2+(ubu1-ub1)**2+(ubu2-ub2)**2])/3
    #p_ch = np.array([float(ub_m)-oub, float(ub_m)-float(ubl_l),
    #                 float(ubu_u)-float(ub_m)])*100/un

    #return oub, float(ub_m), float(ubl_l), float(ubu_u), p_ch


def fig4_multi_final(tmean, conversion, rcpplume, a, b, c,
                tmean1, a1, b1, c1,
                tmean2, a2, b2, c2,
                title=''):
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
        f1 = 3.71 * np.log(ppm1/278) / np.log(2)
        f2 = 3.71 * np.log(ppm2/278) / np.log(2)
        dfdc = (f2-f1)/20
        return dfdc
    fs = 12
    fig, ax = plt.subplots(3, 1, figsize=(0.4375*13, 13*1.5),
                           facecolor='w', edgecolor='k', linewidth=2)
    #fig.delaxes(ax[1, 1])
    # co2 = np.array([390.4, 423.4, 395.8, 550.0])
    '''
    tmean = np.zeros((4))
    tmean[0] = np.mean(tmeans['batch_518'])
    tmean[1] = np.mean(tmeans['batch_520'])
    tmean[2] = np.mean(tmeans['batch_521'])
    tmean[3] = np.mean(tmeans['batch_522'])

    tmean1 = np.zeros((3))
    tmean1[0] = np.mean(tmeans1['present'])
    tmean1[1] = np.mean(tmeans1['plus15'])
    tmean1[2] = np.mean(tmeans1['plus15_lower'])
    tmean1[3] = np.mean(tmeans1['plus15_higher'])

    tmean2 = np.zeros((4))
    tmean2[0] = np.mean(tmeans2['batch_518'])
    tmean2[1] = np.mean(tmeans2['batch_520'])+.05
    tmean2[2] = np.mean(tmeans2['batch_521'])
    tmean2[3] = np.mean(tmeans2['batch_522'])
    '''
    def sigt(x, cv, alpha, beta):
        sigT = np.sqrt((x/beta)**2*cv[1, 1] +
                       (x*alpha/beta**2)**2*cv[0, 0] +
                       2*(x**2*alpha*cv[1, 0]/beta**3))
        return sigT
    
    ensc3 = '#559E54'
    ensc2 = '#183BF0'
    ensc1=  '#9A0794'

    # co2gtc = ppm2gtc(co2)
    rcp_conc = 316.301
    var = ['R95p', 'TX90p', 'WBGT95p']
    let = ['c', 'a', 'b']
    for i in range(3):
        j = 1-i
        if j == -1:
            j = 2
        ax1 = ax[j]
        alpha = a[2-i]
        beta = b[2-i]
        cvs = c[2-i]
        alpha1 = a1[2-i]
        beta1 = b1[2-i]
        cvs1 = c1[2-i]
        alpha2 = a2[2-i]
        beta2 = b2[2-i]
        cvs2 = c2[2-i]
        v = var[2-i]
        l = let[2-i]
        ax1.set_ylim([0, 1.4])
        hfont = {'fontname': 'Arial'}
        ax1.set_ylabel('Global mean temperature anomaly relative to 2006-2015 \
($^{\circ}$C)', fontsize=fs, **hfont)
        ax1.set_xlabel('Cumulative total anthropogenic CO$\mathregular{_{2}}$ \
emissions from 2011 (Gt C)', fontsize=fs, **hfont)
        ax1.set_xticks(np.arange(0, 1100, 100))
        ax1.set_xlim([0, 650])
        # ax1.plot(conversion[:, 2], 0.01951515*(conversion[:, 0]-2010))
        # ax1.plot(conversion[156:162, 2], test1[4:])
        x11 = rcpplume[0, 20]
        y11 = rcpplume[1, 20]
        x22 = rcpplume[0, 3]
        y22 = rcpplume[1, 3]
        g1 = (rcpplume[1, 19]-rcpplume[1, 20])/(rcpplume[0, 19]-rcpplume[0, 20])
        g2 = (rcpplume[1, 4]-rcpplume[1, 3])/(rcpplume[0, 4]-rcpplume[0, 3])
        ax1.axhline((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3,
                    linestyle='--', color='#091034',
                    linewidth=2)
        # bls = ax1.axhspan(tmean1[1]-tmean1[0], tmean2[1]-tmean2[0],
        #                   color='#091034', alpha=0.25)
        ax1.fill_between(rcpplume[0, :], rcpplume[1, :], facecolor='#F7E3D4',
                         edgecolor='#F7E3D4', alpha=1)
        x = np.arange(0, 1010, .01)
    
        # first ensemble
        y = (-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0])
        cns, = ax1.plot(x, y, linewidth=2, color=ensc1)
        y1 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) -
              1.6445*sigt(gtc2df(x), cvs, alpha, beta))
        y2 = ((-alpha*dfdc(x)/beta * ((x)-rcp_conc)+tmean[1]-tmean[0]) +
              1.6445*sigt(gtc2df(x), cvs, alpha, beta))
        ax1.fill_between(x, y1, y2, facecolor=ensc1, edgecolor=ensc1,
                         alpha=.25)
        oub = ((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3 -
               y22) / g2 + x22
        lb = ((tmean[1]-tmean[0]+tmean1[1]-tmean1[0]+tmean2[1]-tmean2[0])/3 -
              y11) / g1 + x11
        ub = x[np.argwhere(np.diff(np.sign(y -
                                   g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu = x[np.argwhere(np.diff(np.sign(y2 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl = x[np.argwhere(np.diff(np.sign(y1 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    
        ax1.axvline(oub, linestyle='--', color='#FD3D12', linewidth=2)
        ax1.axvline(lb, linestyle='--', color='#FD3D12', linewidth=2)
    
        mns = ax1.scatter(rcp_conc, tmean[1]-tmean[0], color=ensc1,
                          marker='o', s=100, lw=3, label='RCP2.6 2090s')
        # second ensemble
        y_1 = (-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0])
        cns1, = ax1.plot(x, y_1, linewidth=2, color=ensc2)
        y1_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) -
                1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
        y2_1 = ((-alpha1*dfdc(x)/beta1 * ((x)-rcp_conc)+tmean1[1]-tmean1[0]) +
                1.6445*sigt(gtc2df(x), cvs1, alpha1, beta1))
        ax1.fill_between(x, y1_1, y2_1, facecolor=ensc2, edgecolor=ensc2,
                         alpha=.25)
        # oub1 = (tmean1[1]-tmean1[0] - y22) / g2 + x22
        # lb1 = (tmean1[1]-tmean1[0] - y11) / g1 + x11
        ub1 = x[np.argwhere(np.diff(np.sign(y_1 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu1 = x[np.argwhere(np.diff(np.sign(y2_1 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl1 = x[np.argwhere(np.diff(np.sign(y1_1 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
    
        mns1 = ax1.scatter(rcp_conc, tmean1[1]-tmean1[0], color=ensc2,
                           marker='o', s=100, lw=3, label='RCP2.6 2090s')
        # third ensemble
        y_2 = (-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0])
        cns2, = ax1.plot(x, y_2, linewidth=2, color=ensc3)
        y1_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) -
                1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
        y2_2 = ((-alpha2*dfdc(x)/beta2 * ((x)-rcp_conc)+tmean2[1]-tmean2[0]) +
                1.6445*sigt(gtc2df(x), cvs2, alpha2, beta2))
        ax1.fill_between(x, y1_2, y2_2, facecolor=ensc3, edgecolor=ensc3,
                         alpha=.25)
        # oub2 = (tmean2[1]-tmean2[0] - y22) / g2 + x22
        # lb2 = (tmean2[1]-tmean2[0] - y11) / g1 + x11
        ub2 = x[np.argwhere(np.diff(np.sign(y_2 -
                                    g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubu2 = x[np.argwhere(np.diff(np.sign(y2_2 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        ubl2 = x[np.argwhere(np.diff(np.sign(y1_2 -
                                     g2*(x-x22)-y22)) != 0).reshape(-1) + 0]
        spr = ax1.axvspan(np.min([ubl, ubl1, ubl2]), np.max([ubu, ubu1, ubu2]),
                          color='#FBAA09', alpha=0.25)
        #spr = ax1.axvspan((ub+ub1+ub2)/3-np.sqrt([(ubl-ub)**2+(ubl1-ub1)**2+(ubl2-ub2)**2])/3, (ub+ub1+ub2)/3+np.sqrt([(ubu-ub)**2+(ubu1-ub1)**2+(ubu2-ub2)**2])/3,
         #                 color='#FBAA09', alpha=0.25)
        #spr = ax1.axvspan((ub+ub1+ub2)/3-1.6445*np.std([ub,ub1,ub2]), (ub+ub1+ub2)/3+1.6445*np.std([ub,ub1,ub2]),
         #                 color='#FBAA09', alpha=0.25)
        # rls = ax1.axvspan(oub2, oub1, color='#FD3D12', alpha=0.25)
        # ax1.axvspan(lb2, lb1, color='#FD3D12', alpha=0.25)
        ubp = ax1.axvline((ub+ub1+ub2)/3, linestyle='--', color='#FBAA09',
                          linewidth=2)
        mns2 = ax1.scatter(rcp_conc, tmean2[1]-tmean2[0], color=ensc3,
                           marker='o', s=100, lw=3, label='RCP2.6 2090s')
    
        # as normal
        ax1.set_yticklabels(ax1.get_yticks(), fontsize=fs, **hfont)
        ax1.set_xticklabels(ax1.get_xticks(), fontsize=fs, **hfont)
        ax1.xaxis.set_tick_params(width=2, length=7.5)
        ax1.yaxis.set_tick_params(width=2, length=7.5)
        # ax2.xaxis.set_tick_params(width=2, length=7.5)
        ax1.axvline(0, linewidth=2.5, color='k')
        ax1.axvline(1000, linewidth=2.5, color='k')
        ax1.axhline(0, linewidth=2.5, color='k')
        ax1.axhline(1.4, linewidth=2.5, color='k')
    
        pink_patch = mpatches.Patch(color='#F7E3D4', label='RCP range')
        red_line = mlines.Line2D([], [], color='#FD3D12', linestyle='--',
                                 linewidth=2, label='Carbon budget bound')
        black_line = mlines.Line2D([], [], color='#091034', linestyle='--',
                                   linewidth=2, label='1.5$^{\circ}$C warming')
        b_spr = mpatches.Patch(color=ensc1, alpha=0.25, linewidth=0)
        b_spr1 = mpatches.Patch(color=ensc2, alpha=0.25, linewidth=0)
        b_spr2 = mpatches.Patch(color=ensc3, alpha=0.25, linewidth=0)
        ax1.legend(handles=[pink_patch, black_line, mns, (cns, b_spr), mns1,
                            (cns1, b_spr1), mns2, (cns2, b_spr2),
                            red_line, (ubp, spr)],
                   labels=['RCP range', '1.5$^{\circ}$C warming',
                           'RCP2.6 2090s HadAM3P', 'Constant '+v+' HadAM3P',
                           'RCP2.6 2090s MIROC5', 'Constant '+v+' MIROC5',
                           'RCP2.6 2090s CAM4', 'Constant '+v+' CAM4',
                           'Carbon budget bound', 'New upper bound'],
                   handlelength=3,
                   loc=1, scatterpoints=1, prop={'family': 'Arial', 'size':fs-4},
                   bbox_to_anchor=(.995, .995), frameon=True, framealpha=1)
        ax1.set_title(l, loc='left', fontname='Arial', fontsize=fs+2,
                      fontweight='bold', x=-.15, y=1)
    plt.savefig('/home/bakerh/Downloads/fig4.png', dpi=300, bbox=np.array([[0,0],[3.4647, 11.8785]]))
    #plt.title(title, y=1.01, fontsize=20, **hfont)
    # plt.subplots_adjust(hspace=0, wspace=0.05, top=.95, bottom=0.1, left=.05,
    #                   right=.97)
    un = oub - lb
    ub_m = (ub+ub1+ub2)/3
    ubl_l = np.min([ubl, ubl1, ubl2])
    ubu_u = np.max([ubu, ubu1, ubu2])
    #ubl_l = (ub+ub1+ub2)/3-np.sqrt([(ubl-ub)**2+(ubl1-ub1)**2+(ubl2-ub2)**2])/3
    #ubu_u = (ub+ub1+ub2)/3+np.sqrt([(ubu-ub)**2+(ubu1-ub1)**2+(ubu2-ub2)**2])/3
    p_ch = np.array([float(ub_m)-oub, float(ub_m)-float(ubl_l),
                     float(ubu_u)-float(ub_m)])*100/un

    return oub, float(ub_m), float(ubl_l), float(ubu_u), p_ch


def boxplot_multi(globprcnts, prcnts, globprcntsm, prcntsm):
    from matplotlib import pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    hfont = {'fontname': 'Arial'}
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([0.5, 1, 1.5, 2, 2.5, 3,
                  4.5, 5, 5.5, 6, 6.5, 7,
                  8.5, 9, 9.5, 10, 10.5, 11,
                  14, 14.5, 15, 15.5, 16, 16.5,
                  16, 16.5, 17, 17.5, 18, 18.5,
                  21.5, 22, 22.5, 23, 23.5, 24,
                  25.5, 26, 26.5, 27, 27.5, 28,
                  31, 31.5, 32, 32.5, 33, 33.5,
                  35, 35.5, 36, 36.5, 37, 37.5,
                  39, 39.5, 40, 40.5, 41, 41.5,
                  44.5, 45, 45.5, 46, 46.5, 47,
                  48.5, 49, 49.5, 50, 50.5, 51])
    x[24:] += 2
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

    globprcnts2m = np.copy(globprcntsm)
    globprcnts2m[0, 1:, 0] -= globprcntsm[0, 0, 0]
    globprcnts2m = np.vstack((globprcnts2m[0, 2], globprcnts2m[0, 1],
                             globprcnts2m[0, 3]))
    globprcnts3m = np.copy(globprcntsm)
    globprcnts3m[2, 1:, 0] -= globprcntsm[2, 0, 0]
    globprcnts3m = np.vstack((globprcnts3m[2, 2], globprcnts3m[2, 1],
                             globprcnts3m[2, 3]))
    prcnts2m = np.copy(prcntsm)
    prcnts2m[1:, :, 0] -= prcntsm[0, :, 0]
    prcnts2m = prcnts2m[1:, :]
    holdm = np.copy(prcnts2m[0, :])
    prcnts2m[0, :] = prcnts2m[1, :]
    prcnts2m[1, :] = holdm

    labels = ['', '', '', '$\mathregular{\overline{T}_{GLOB}}$', '', '', '', '', '', '$\mathregular{\overline{T}_{NH\_ET}}$', '', '', '', '', '','$\mathregular{\overline{T}_{TROP}}$', '', '', '', '', '',
              '$\mathregular{TX90p_{NH\_ET}}$', '', '', '', '', '', '$\mathregular{TX90p_{TROP}}$', '', '', '', '', '', '$\mathregular{WBGT95p_{NH\_ET}}$', '', '','', '', '',
              '$\mathregular{WBGT95p_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{R}_{GLOB}}$', '', '', '', '', '', '$\mathregular{\overline{R}_{NH\_ET}}$', '', '', '', '', '',
              '$\mathregular{\overline{R}_{TROP}}$', '', '', '', '', '', '$\mathregular{R95p_{NH\_ET}}$', '', '', '', '', '', '$\mathregular{R95p_{TROP}}$', '', '']

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

    ym = np.concatenate((globprcnts2m[:, 0].flatten(order='F'),
                        prcnts2m[:, :6, 0].flatten(order='F')))
    yerrsm = np.array([np.concatenate((globprcnts2m[:, 1].flatten(order='F'),
                                      prcnts2m[:, :6, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2m[:, 2].flatten(order='F'),
                                      prcnts2m[:, :6, 2].flatten(order='F')))])
    y1m = np.concatenate((globprcnts3m[:, 0].flatten(order='F'),
                        prcnts2m[:, 6:, 0].flatten(order='F')))
    yerrs1m = np.array([np.concatenate((globprcnts3m[:, 1].flatten(order='F'),
                                      prcnts2m[:, 6:, 1].flatten(order='F'))),
                      np.concatenate((globprcnts3m[:, 2].flatten(order='F'),
                                      prcnts2m[:, 6:, 2].flatten(order='F')))])
    color = ['b', 'k', 'r']
    pp = []
    for i in range(3):
        p = ax1.errorbar(x[i:18:6], y[i:9:3], yerr=yerrs[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(p[2])
    ax2 = ax1.twinx()
    for i in range(3):
        ax2.errorbar(x[18+i:30:6], y[9+i:15:3], yerr=yerrs[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    for i in range(3):
        ax3.errorbar(x[30+i:42:6], y[15+i::3], yerr=yerrs[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax4 = ax1.twinx()
    for i in range(3):
        ax4.errorbar(x[42+i:60:6], y1[i:9:3]*30, yerr=yerrs1[:, i:9:3]*30,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax5 = ax1.twinx()
    for i in range(3):
        ax5.errorbar(x[60+i::6], y1[9+i::3], yerr=yerrs1[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)

    for i in range(3):
        pm = ax1.errorbar(x[i+3:18:6], ym[i:9:3], yerr=yerrsm[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(pm[2])
        pm[-1][0].set_linestyle('--')

    for i in range(3):
        pn=ax2.errorbar(x[18+i+3:30:6], ym[9+i:15:3], yerr=yerrsm[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pn[-1][0].set_linestyle('--')

    for i in range(3):
        po=ax3.errorbar(x[30+i+3:42:6], ym[15+i::3], yerr=yerrsm[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        po[-1][0].set_linestyle('--')

    for i in range(3):
        pq=ax4.errorbar(x[42+i+3:60:6], y1m[i:9:3]*30, yerr=yerrs1m[:, i:9:3]*30,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pq[-1][0].set_linestyle('--')

    for i in range(3):
        pr=ax5.errorbar(x[60+i+3::6], y1m[9+i::3], yerr=yerrs1m[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pr[-1][0].set_linestyle('--')

    ax2.spines['left'].set_position(('data', 13.5))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 23))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax4.spines['left'].set_position(('data', 32.5))
    ax4.spines['right'].set_color('none')
    ax4.spines['bottom'].set_position(('axes', 0.0))
    ax4.spines['top'].set_color('none')
    ax4.spines['bottom'].set_smart_bounds(True)
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')

    ax5.spines['left'].set_position(('data', 46))
    ax5.spines['right'].set_color('none')
    ax5.spines['bottom'].set_position(('axes', 0.0))
    ax5.spines['top'].set_color('none')
    ax5.spines['bottom'].set_smart_bounds(True)
    ax5.xaxis.set_ticks_position('bottom')
    ax5.yaxis.set_ticks_position('left')
    ax5.yaxis.set_label_position('left')

    ax1.set_xlim([0, 53.5])
    ax1.set_ylim([.6, 1.8])
    ax2.set_ylim([3.75, 17.75])
    ax3.set_ylim([.175, 2.2])
    ax4.set_ylim([-1.24, 3.75])
    ax5.set_ylim([-0.025, .475])
    ax1.set_ylabel('Temperature anomaly ($^{\circ}$C)', fontsize=20, **hfont)
    ax2.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)', fontsize=20,
                   **hfont)
    ax3.set_ylabel('Anomaly ($^{\circ}$C months)',
                   fontsize=20, labelpad=-5, **hfont)
    ax4.set_ylabel('Precipitation anomaly (mm month$\mathregular{^{-1}}$)',
                   fontsize=20, labelpad=-5, **hfont)
    ax5.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)',
                   fontsize=20, labelpad=-5, **hfont)

    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation='vertical', **hfont)
    ax2.set_yticklabels(ax2.get_yticks(), fontsize=20, **hfont)
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=20, **hfont)
    ax4.set_yticklabels(ax4.get_yticks(), fontsize=20, **hfont)
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
    ax5.set_yticklabels(ax5.get_yticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax2.yaxis.set_tick_params(width=2, length=7.5)
    ax3.yaxis.set_tick_params(width=2, length=7.5)
    ax4.yaxis.set_tick_params(width=2, length=7.5)
    ax5.yaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(13.5, linewidth=1.5, color='k')
    ax1.axvline(23, linewidth=1.5, color='k')
    ax1.axvline(32.5, linewidth=1.5, color='k')
    ax1.axvline(46, linewidth=1.5, color='k')
    ax1.axvline(53.5, linewidth=2.5, color='k')
    ax1.axhline(0.6, linewidth=2.5, color='k')
    ax1.axhline(1.8, linewidth=2.5, color='k')
    ax1.legend(handles=pp, labels=['HadAM3P Low',
                                   'HadAM3P B.E.',
                                   'HadAM3P High',
                                   'MIROC5 Low',
                                   'MIROC5 B.E.',
                                   'MIROC5 High'],
               loc=2, scatterpoints=3, prop={'family': 'Arial'},
               bbox_to_anchor=(.0025, .89075, 1., .102))
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.2, left=.05,
                        right=.99)


def boxplot_multi_forcing(globprcnts, prcnts, globprcntsm, prcntsm):
    from matplotlib import pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    hfont = {'fontname': 'Arial'}
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([0.5, 1, 1.5, 2, 2.5, 3,
                  4.5, 5, 5.5, 6, 6.5, 7,
                  8.5, 9, 9.5, 10, 10.5, 11,
                  14, 14.5, 15, 15.5, 16, 16.5,
                  16, 16.5, 17, 17.5, 18, 18.5,
                  21.5, 22, 22.5, 23, 23.5, 24,
                  25.5, 26, 26.5, 27, 27.5, 28,
                  31, 31.5, 32, 32.5, 33, 33.5,
                  35, 35.5, 36, 36.5, 37, 37.5,
                  39, 39.5, 40, 40.5, 41, 41.5,
                  44.5, 45, 45.5, 46, 46.5, 47,
                  48.5, 49, 49.5, 50, 50.5, 51])
    x[24:] += 2
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

    globprcnts2m = np.copy(globprcntsm)
    globprcnts2m[0, 1:, 0] -= globprcntsm[0, 0, 0]
    globprcnts2m = np.vstack((globprcnts2m[0, 2], globprcnts2m[0, 1],
                             globprcnts2m[0, 3]))
    globprcnts3m = np.copy(globprcntsm)
    globprcnts3m[2, 1:, 0] -= globprcntsm[2, 0, 0]
    globprcnts3m = np.vstack((globprcnts3m[2, 2], globprcnts3m[2, 1],
                             globprcnts3m[2, 3]))
    prcnts2m = np.copy(prcntsm)
    prcnts2m[1:, :, 0] -= prcntsm[0, :, 0]
    prcnts2m = prcnts2m[1:, :]
    holdm = np.copy(prcnts2m[0, :])
    prcnts2m[0, :] = prcnts2m[1, :]
    prcnts2m[1, :] = holdm

    labels = ['', '', '', '$\mathregular{\overline{T}_{GLOB}}$', '', '', '', '', '', '$\mathregular{\overline{T}_{NH\_ET}}$', '', '', '', '', '','$\mathregular{\overline{T}_{TROP}}$', '', '', '', '', '',
              '$\mathregular{\overline{TX90p}_{NH\_ET}}$', '', '', '', '', '', '$\mathregular{\overline{TX90p}_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{WBGT95p}_{NH\_ET}}$', '', '','', '', '',
              '$\mathregular{\overline{WBGT95p}_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{R}_{GLOB}}$', '', '', '', '', '', '$\mathregular{\overline{R}_{NH\_ET}}$', '', '', '', '', '',
              '$\mathregular{\overline{R}_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{R95p}_{NH\_ET}}$', '', '', '', '', '', '$\mathregular{\overline{R95p}_{TROP}}$', '', '']

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

    ym = np.concatenate((globprcnts2m[:, 0].flatten(order='F'),
                        prcnts2m[:, :6, 0].flatten(order='F')))
    yerrsm = np.array([np.concatenate((globprcnts2m[:, 1].flatten(order='F'),
                                      prcnts2m[:, :6, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2m[:, 2].flatten(order='F'),
                                      prcnts2m[:, :6, 2].flatten(order='F')))])
    y1m = np.concatenate((globprcnts3m[:, 0].flatten(order='F'),
                        prcnts2m[:, 6:, 0].flatten(order='F')))
    yerrs1m = np.array([np.concatenate((globprcnts3m[:, 1].flatten(order='F'),
                                      prcnts2m[:, 6:, 1].flatten(order='F'))),
                      np.concatenate((globprcnts3m[:, 2].flatten(order='F'),
                                      prcnts2m[:, 6:, 2].flatten(order='F')))])
    color = ['b', 'k', 'r']
    pp = []
    for i in range(3):
        p = ax1.errorbar(x[i:18:6], y[i:9:3], yerr=yerrs[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(p[2])
    ax2 = ax1.twinx()
    for i in range(3):
        ax2.errorbar(x[18+i:30:6], y[9+i:15:3], yerr=yerrs[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    for i in range(3):
        ax3.errorbar(x[30+i:42:6], y[15+i::3], yerr=yerrs[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax4 = ax1.twinx()
    for i in range(3):
        ax4.errorbar(x[42+i:60:6], y1[i:9:3]*30, yerr=yerrs1[:, i:9:3]*30,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax5 = ax1.twinx()
    for i in range(3):
        ax5.errorbar(x[60+i::6], y1[9+i::3], yerr=yerrs1[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)

    for i in range(3):
        pm = ax1.errorbar(x[i+3:18:6], ym[i:9:3], yerr=yerrsm[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(pm[2])
        pm[-1][0].set_linestyle('--')

    for i in range(3):
        pn=ax2.errorbar(x[18+i+3:30:6], ym[9+i:15:3], yerr=yerrsm[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pn[-1][0].set_linestyle('--')

    for i in range(3):
        po=ax3.errorbar(x[30+i+3:42:6], ym[15+i::3], yerr=yerrsm[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        po[-1][0].set_linestyle('--')

    for i in range(3):
        pq=ax4.errorbar(x[42+i+3:60:6], y1m[i:9:3]*30, yerr=yerrs1m[:, i:9:3]*30,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pq[-1][0].set_linestyle('--')

    for i in range(3):
        pr=ax5.errorbar(x[60+i+3::6], y1m[9+i::3], yerr=yerrs1m[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pr[-1][0].set_linestyle('--')

    ax2.spines['left'].set_position(('data', 13.5))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 23))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax4.spines['left'].set_position(('data', 32.5))
    ax4.spines['right'].set_color('none')
    ax4.spines['bottom'].set_position(('axes', 0.0))
    ax4.spines['top'].set_color('none')
    ax4.spines['bottom'].set_smart_bounds(True)
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')

    ax5.spines['left'].set_position(('data', 46))
    ax5.spines['right'].set_color('none')
    ax5.spines['bottom'].set_position(('axes', 0.0))
    ax5.spines['top'].set_color('none')
    ax5.spines['bottom'].set_smart_bounds(True)
    ax5.xaxis.set_ticks_position('bottom')
    ax5.yaxis.set_ticks_position('left')
    ax5.yaxis.set_label_position('left')

    ax1.set_xlim([0, 53.5])
    ax1.set_ylim([-.1, 0.8])
    ax2.set_ylim([-.8, 5.2])
    ax3.set_ylim([-.1, .95])
    ax4.set_ylim([-2.5, 2.75])
    ax5.set_ylim([-0.12, .29])
    ax1.set_ylabel('Temperature anomaly ($^{\circ}$C)', fontsize=20, **hfont)
    ax2.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)', fontsize=20,
                   **hfont)
    ax3.set_ylabel('Anomaly ($^{\circ}$C months)',
                   fontsize=20, labelpad=-5, **hfont)
    ax4.set_ylabel('Precipitation anomaly (mm month$\mathregular{^{-1}}$)',
                   fontsize=20, labelpad=-5, **hfont)
    ax5.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)',
                   fontsize=20, labelpad=-5, **hfont)

    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation='vertical', **hfont)
    ax2.set_yticklabels(ax2.get_yticks(), fontsize=20, **hfont)
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=20, **hfont)
    ax4.set_yticklabels(ax4.get_yticks(), fontsize=20, **hfont)
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
    ax5.set_yticklabels(ax5.get_yticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax2.yaxis.set_tick_params(width=2, length=7.5)
    ax3.yaxis.set_tick_params(width=2, length=7.5)
    ax4.yaxis.set_tick_params(width=2, length=7.5)
    ax5.yaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(13.5, linewidth=1.5, color='k')
    ax1.axvline(23, linewidth=1.5, color='k')
    ax1.axvline(32.5, linewidth=1.5, color='k')
    ax1.axvline(46, linewidth=1.5, color='k')
    ax1.axvline(53.5, linewidth=2.5, color='k')
    ax1.axhline(-.1, linewidth=2.5, color='k')
    ax1.axhline(0.8, linewidth=2.5, color='k')
    ax1.legend(handles=pp, labels=['HadAM3P Low',
                                   'HadAM3P B.E.',
                                   'HadAM3P High',
                                   'MIROC5 Low',
                                   'MIROC5 B.E.',
                                   'MIROC5 High'],
               loc=2, scatterpoints=3, prop={'family': 'Arial'},
               bbox_to_anchor=(.0025, .89075, 1., .102))
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.2, left=.05,
                        right=.99)


def boxplot_multi_forcing2(globprcnts, prcnts, globprcntsm, prcntsm):
    from matplotlib import pyplot as plt
    from matplotlib.ticker import FormatStrFormatter
    hfont = {'fontname': 'Arial'}
    fig, ax1 = plt.subplots(facecolor='w', edgecolor='k', linewidth=2)
    x = np.array([0.5, 1, 1.5, 2, 2.5, 3,
                  4.5, 5, 5.5, 6, 6.5, 7,
                  8.5, 9, 9.5, 10, 10.5, 11,
                  14, 14.5, 15, 15.5, 16, 16.5,
                  16, 16.5, 17, 17.5, 18, 18.5,
                  21.5, 22, 22.5, 23, 23.5, 24,
                  25.5, 26, 26.5, 27, 27.5, 28,
                  31, 31.5, 32, 32.5, 33, 33.5,
                  35, 35.5, 36, 36.5, 37, 37.5,
                  39, 39.5, 40, 40.5, 41, 41.5,
                  44.5, 45, 45.5, 46, 46.5, 47,
                  48.5, 49, 49.5, 50, 50.5, 51])
    x[24:] += 2
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

    globprcnts2m = np.copy(globprcntsm)
    globprcnts2m[0, 1:, 0] -= globprcntsm[0, 0, 0]
    globprcnts2m = np.vstack((globprcnts2m[0, 2], globprcnts2m[0, 1],
                             globprcnts2m[0, 3]))
    globprcnts3m = np.copy(globprcntsm)
    globprcnts3m[2, 1:, 0] -= globprcntsm[2, 0, 0]
    globprcnts3m = np.vstack((globprcnts3m[2, 2], globprcnts3m[2, 1],
                             globprcnts3m[2, 3]))
    prcnts2m = np.copy(prcntsm)
    prcnts2m[1:, :, 0] -= prcntsm[0, :, 0]
    prcnts2m = prcnts2m[1:, :]
    holdm = np.copy(prcnts2m[0, :])
    prcnts2m[0, :] = prcnts2m[1, :]
    prcnts2m[1, :] = holdm

    labels = ['', '', '', '$\mathregular{\overline{T}_{GLOB}}$', '', '', '', '', '', '$\mathregular{\overline{T}_{NH\_ET}}$', '', '', '', '', '','$\mathregular{\overline{T}_{TROP}}$', '', '', '', '', '',
              '$\mathregular{\overline{TX90p}_{NH\_ET}}$', '', '', '', '', '', '$\mathregular{\overline{TX90p}_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{WBGT95p}_{NH\_ET}}$', '', '','', '', '',
              '$\mathregular{\overline{WBGT95p}_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{R}_{GLOB}}$', '', '', '', '', '', '$\mathregular{\overline{R}_{NH\_ET}}$', '', '', '', '', '',
              '$\mathregular{\overline{R}_{TROP}}$', '', '', '', '', '', '$\mathregular{\overline{R95p}_{NH\_ET}}$', '', '', '', '', '', '$\mathregular{\overline{R95p}_{TROP}}$', '', '']

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

    ym = np.concatenate((globprcnts2m[:, 0].flatten(order='F'),
                        prcnts2m[:, :6, 0].flatten(order='F')))
    yerrsm = np.array([np.concatenate((globprcnts2m[:, 1].flatten(order='F'),
                                      prcnts2m[:, :6, 1].flatten(order='F'))),
                      np.concatenate((globprcnts2m[:, 2].flatten(order='F'),
                                      prcnts2m[:, :6, 2].flatten(order='F')))])
    y1m = np.concatenate((globprcnts3m[:, 0].flatten(order='F'),
                        prcnts2m[:, 6:, 0].flatten(order='F')))
    yerrs1m = np.array([np.concatenate((globprcnts3m[:, 1].flatten(order='F'),
                                      prcnts2m[:, 6:, 1].flatten(order='F'))),
                      np.concatenate((globprcnts3m[:, 2].flatten(order='F'),
                                      prcnts2m[:, 6:, 2].flatten(order='F')))])
    color = ['b', 'k', 'r']
    pp = []
    for i in range(3):
        p = ax1.errorbar(x[i:18:6], y[i:9:3], yerr=yerrs[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(p[2])
    ax2 = ax1.twinx()
    for i in range(3):
        ax2.errorbar(x[18+i:30:6], y[9+i:15:3], yerr=yerrs[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax3 = ax1.twinx()
    for i in range(3):
        ax3.errorbar(x[30+i:42:6], y[15+i::3], yerr=yerrs[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax4 = ax1.twinx()
    for i in range(3):
        ax4.errorbar(x[42+i:60:6], y1[i:9:3]*30, yerr=yerrs1[:, i:9:3]*30,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
    ax5 = ax1.twinx()
    for i in range(3):
        ax5.errorbar(x[60+i::6], y1[9+i::3], yerr=yerrs1[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)

    for i in range(3):
        pm = ax1.errorbar(x[i+3:18:6], ym[i:9:3], yerr=yerrsm[:, i:9:3],
                         color=color[i], fmt='x', elinewidth=2, capthick=2,
                         capsize=5, ms=10, mew=2)
        pp.append(pm[2])
        pm[-1][0].set_linestyle('--')

    for i in range(3):
        pn=ax2.errorbar(x[18+i+3:30:6], ym[9+i:15:3], yerr=yerrsm[:, 9+i:15:3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pn[-1][0].set_linestyle('--')

    for i in range(3):
        po=ax3.errorbar(x[30+i+3:42:6], ym[15+i::3], yerr=yerrsm[:, 15+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        po[-1][0].set_linestyle('--')

    for i in range(3):
        pq=ax4.errorbar(x[42+i+3:60:6], y1m[i:9:3]*30, yerr=yerrs1m[:, i:9:3]*30,
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pq[-1][0].set_linestyle('--')

    for i in range(3):
        pr=ax5.errorbar(x[60+i+3::6], y1m[9+i::3], yerr=yerrs1m[:, 9+i::3],
                     color=color[i], fmt='x', elinewidth=2, capthick=2,
                     capsize=5, ms=10, mew=2)
        pr[-1][0].set_linestyle('--')

    ax2.spines['left'].set_position(('data', 13.5))
    ax2.spines['right'].set_color('none')
    ax2.spines['bottom'].set_position(('axes', 0.0))
    ax2.spines['top'].set_color('none')
    ax2.spines['bottom'].set_smart_bounds(True)
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.yaxis.set_label_position('left')

    ax3.spines['left'].set_position(('data', 23))
    ax3.spines['right'].set_color('none')
    ax3.spines['bottom'].set_position(('axes', 0.0))
    ax3.spines['top'].set_color('none')
    ax3.spines['bottom'].set_smart_bounds(True)
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.yaxis.set_label_position('left')

    ax4.spines['left'].set_position(('data', 32.5))
    ax4.spines['right'].set_color('none')
    ax4.spines['bottom'].set_position(('axes', 0.0))
    ax4.spines['top'].set_color('none')
    ax4.spines['bottom'].set_smart_bounds(True)
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.yaxis.set_label_position('left')

    ax5.spines['left'].set_position(('data', 46))
    ax5.spines['right'].set_color('none')
    ax5.spines['bottom'].set_position(('axes', 0.0))
    ax5.spines['top'].set_color('none')
    ax5.spines['bottom'].set_smart_bounds(True)
    ax5.xaxis.set_ticks_position('bottom')
    ax5.yaxis.set_ticks_position('left')
    ax5.yaxis.set_label_position('left')

    ax1.set_xlim([0, 53.5])
    ax1.set_ylim([.6, 1.625])
    ax2.set_ylim([3.75, 15.25])
    ax3.set_ylim([.175, 1.975])
    ax4.set_ylim([-1.3, 3.5])
    ax5.set_ylim([-0.025, .45])
    ax1.set_ylabel('Temperature anomaly ($^{\circ}$C)', fontsize=20, **hfont)
    ax2.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)', fontsize=20,
                   **hfont)
    ax3.set_ylabel('Anomaly ($^{\circ}$C months)',
                   fontsize=20, labelpad=-5, **hfont)
    ax4.set_ylabel('Precipitation anomaly (mm month$\mathregular{^{-1}}$)',
                   fontsize=20, labelpad=-5, **hfont)
    ax5.set_ylabel('Anomaly (days season$\mathregular{^{-1}}$)',
                   fontsize=20, labelpad=-5, **hfont)

    ax1.set_xticks(x)
    ax1.set_yticklabels(ax1.get_yticks(), fontsize=20, **hfont)
    ax1.set_xticklabels(labels, fontsize=20,
                        rotation='vertical', **hfont)
    ax2.set_yticklabels(ax2.get_yticks(), fontsize=20, **hfont)
    ax3.set_yticklabels(ax3.get_yticks(), fontsize=20, **hfont)
    ax4.set_yticklabels(ax4.get_yticks(), fontsize=20, **hfont)
    ax4.yaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
    ax5.set_yticklabels(ax5.get_yticks(), fontsize=20, **hfont)
    ax1.xaxis.set_tick_params(width=2, length=7.5)
    ax1.yaxis.set_tick_params(width=2, length=7.5)
    ax2.yaxis.set_tick_params(width=2, length=7.5)
    ax3.yaxis.set_tick_params(width=2, length=7.5)
    ax4.yaxis.set_tick_params(width=2, length=7.5)
    ax5.yaxis.set_tick_params(width=2, length=7.5)
    ax1.axvline(0, linewidth=2.5, color='k')
    ax1.axvline(13.5, linewidth=1.5, color='k')
    ax1.axvline(23, linewidth=1.5, color='k')
    ax1.axvline(32.5, linewidth=1.5, color='k')
    ax1.axvline(46, linewidth=1.5, color='k')
    ax1.axvline(53.5, linewidth=2.5, color='k')
    ax1.axhline(0.6, linewidth=2.5, color='k')
    ax1.axhline(1.8, linewidth=2.5, color='k')
    ax1.legend(handles=pp, labels=['HadAM3P Low',
                                   'HadAM3P B.E.',
                                   'HadAM3P High',
                                   'MIROC5 Low',
                                   'MIROC5 B.E.',
                                   'MIROC5 High'],
               loc=2, scatterpoints=3, prop={'family': 'Arial'},
               bbox_to_anchor=(.0025, .89075, 1., .102))
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.2, left=.05,
                        right=.99)
