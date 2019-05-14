#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 10 15:09:45 2017

@author: bakerh
"""

import numpy as np


def hadedge(potemp, z, temp, strm):
    def lapserate(t, z):
        '''
        dT = np.zeros((37, 64))
        dz = np.zeros((37, 64))
        for i in range(36):
            dT[i, :] = t[i+1, :] - t[i, :]
        for i in range(36):
            dz[i, :] = z[i+1, :] - z[i, :]
        lapse = -1000 * dT[0:-1] / dz[0:-1]
        '''
        lpse = -1000*np.gradient(t, z)
        # zonalplot(lapse, sigma[0:-1], lat, 'Lapse rate')
        return lpse
    hadext = np.zeros((2, 307))
    for i in range(307):
        lat_max = np.remainder(np.argmax(strm[i, :25]), 64)
        lat_min = np.remainder(np.argmin(strm[i, :, 32:]), 32) + 32
        lapse = lapserate(temp[i, :, lat_max], z[i, :, lat_max])
        lapse1 = lapserate(temp[i, :, lat_min], z[i, :, lat_min])
        h = np.interp(2, lapse, z[0, :, lat_max])
        h1 = np.interp(2, lapse1, z[0, :, lat_min])
        dv = potemp[i, -1, lat_max] - np.interp(2, lapse,
                                                potemp[i, :, lat_max])
        dv1 = potemp[i, -1, lat_min] - np.interp(2, lapse1,
                                                 potemp[i, :, lat_min])
        hadext[0, i] = (-1*dv*h)**.25
        hadext[1, i] = (-1*dv1*h1)**.25
    return hadext


def rossby(strm, u, lat):
    w = 7.2921e-5
    r = 6.371e6
    # uy = np.gradient(u, r*lat*np.pi/180, axis=2)
    uy = np.gradient(u*np.cos(lat*np.pi/180), lat*np.pi/180, axis=2)
    uy /= r*np.cos(lat*np.pi/180)
    ro = np.zeros((307, 2))
    for i in range(307):
        lat_max = np.remainder(np.argmax(strm[i, :25]), 64)
        duy = np.mean(uy[i, 18:21, lat_max])
        f = 2 * w * np.sin(lat_max*np.pi/180)
        ro[i, 0] = duy / f
        lat_min = np.remainder(np.argmin(strm[i, :, 32:]), 32) + 32
        duy = np.mean(uy[i, 18:21, lat_min])
        f = 2 * w * np.sin(lat_min*np.pi/180)
        ro[i, 1] = duy / f
    return ro


def zonalhadsens(data1, p, sigma, lat, control, csigma, clat):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    data = np.copy(data1)
    fig, axs = plt.subplots(2, 2, sharey='row', sharex='col', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Hadley Cell Sensitivity", size=24, y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    data[0, :] = data[0, :] / 1e10
    caxismin, caxismax, ctrs = -3, 3, np.array([-3,-2.6,-2.2,-1.8,-1.4,-1,-.6,-.2,.2,.6,1,1.4,1.8,2.2,2.6,3])
    caxismin1, caxismax1, ctr1 = -4, 4, np.linspace(-4, 4, 17)

    control = control / 1e11
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[1, 0, :], ctr1,
                              extend='both',
                              cmap=newcmap, vmin=caxismin1, vmax=caxismax1)
    axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                      colors='k', linewidths=1.5)
    #  axs[0, 0].clabel(plot1, [-1,-.5,.5,1], inline=True, inline_spacing=-3,
    #                   fontsize=10,
    #                   fmt='%.1f')

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, 0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                      colors='k', linewidths=1.5)

    plot = axs[1, 0].contourf(meshlat, meshsigma, data[1, 1, :], ctr1,
                              extend='both', cmap=newcmap, vmin=caxismin1,
                              vmax=caxismax1)
    axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                      colors='k', linewidths=1.5)

    plot2 = axs[1, 1].contourf(meshlat, meshsigma, data[0, 1, :], ctrs,
                               extend='both', cmap=newcmap, vmin=caxismin,
                               vmax=caxismax)
    axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                      colors='k', linewidths=1.5)

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0.1)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.2, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(), fontsize=16)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(), fontsize=16)
            meshsigma2 = np.copy(meshsigma)
            meshsigma1 = np.copy(meshsigma)
            # axs[i, j].contour(meshlat, meshsigma, p[i, 1-j, :], [0.005],
            # linewidths=2)
            for a in range(9):
                for b in range(34):
                    if p[1-j, i, a, b] <= 0.005:
                        meshsigma2[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma2, marker='o', color='gray',
                              s=0.5)
            for a in range(9):
                for b in range(34):
                    if p[1-j, i, a, b] > 0.005:
                        meshsigma1[a, b] = np.nan
            axs[i, j].scatter(meshlat, meshsigma1, marker='o', color='black',
                              s=5)
    plt.subplots_adjust(hspace=0.08, wspace=0.045, top=.95, bottom=0.15,
                        left=.1, right=.9)
    cbar_ax = fig.add_axes([0.1, 0.1, 0.391, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='uniform', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin1, caxismax1, 9))
    b.set_label(label='Poleward Hadley Cell extent shift (deg)', size=20)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontsize=16)
    cbar_ax1 = fig.add_axes([0.509, 0.1, 0.391, 0.015])
    a = fig.colorbar(plot2, cax=cbar_ax1, spacing='uniform', extend='max',
                     orientation='horizontal',
                     ticks=np.array([-3,-2.2,-1.4,-.6,0,.6,1.4,2.2,3]))
    a.set_label(label='Change in strength of Hadley Cell (kg s$\mathregular{^{-1}}$ x 10$\mathregular{^{10}}$)', size=20)
    cl = plt.getp(cbar_ax1, 'xmajorticklabels')
    plt.setp(cl, fontsize=16)
    cols = ['Latitude', 'Magnitude']
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, col in zip(axs[0], cols):
        ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24)
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=24)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    plt.show()


def zonalhaditcz(data1, p, sigma, lat, control, csigma, clat, colormax=5):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
    data = np.copy(data1)
    fig, axs = plt.subplots(2, 1, sharey='row', sharex='col', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Hadley Cell SNR Sensitivity", size=24, y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-colormax, colormax, 17)

    control = control / 1e11
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot = axs[0].contourf(meshlat, meshsigma, data[0, :], ctrs,
                           extend='both',
                           cmap=newcmap, vmin=-colormax, vmax=colormax)
    plot1 = axs[0].contour(cmeshlat, cmeshsigma, control, cctrs,
                           colors='k', linewidths=1.5)

    plot = axs[1].contourf(meshlat, meshsigma, data[1, :], ctrs,
                           extend='both', cmap=newcmap, vmin=-colormax,
                           vmax=colormax)
    plot1 = axs[1].contour(cmeshlat, cmeshsigma, control, cctrs,
                           colors='k', linewidths=1.5)

    for i in range(2):
        axs[i].invert_yaxis()
        axs[i].set_xlim([lat[0], lat[-1]])
        axs[i].set_ylim(1, 0.1)
        axs[i].set_xticks(np.arange(-75, 90, 15))
        axs[i].set_yticks(np.arange(0.2, 1.2, 0.2))
        axs[i].set_xticklabels(axs[i].get_xticks(), fontsize=20)
        axs[i].set_yticklabels(axs[i].get_yticks(), fontsize=20)
        meshsigma2 = np.copy(meshsigma)
        meshsigma1 = np.copy(meshsigma)
        # axs[i, j].contour(meshlat, meshsigma, p[i, 1-j, :], [0.005],
        # linewidths=2)
        for a in range(9):
            for b in range(34):
                if p[i, a, b] <= 0.005:
                    meshsigma2[a, b] = np.nan
        axs[i].scatter(meshlat, meshsigma2, marker='o', color='gray',
                       s=0.5)
        for a in range(9):
            for b in range(34):
                if p[i, a, b] > 0.005:
                    meshsigma1[a, b] = np.nan
        axs[i].scatter(meshlat, meshsigma1, marker='o', color='black',
                       s=5)
    plt.subplots_adjust(hspace=0.079, wspace=0.15, top=.95, bottom=0.05,
                        left=.2, right=.8)
    cbar_ax = fig.add_axes([0.85, 0.2, 0.01, 0.6])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='vertical',
                     ticks=np.linspace(-colormax, colormax, 9))
    #b.set_label(label='Hadley cell to ITCZ distance change (deg)', size=20)
    b.set_label(label="Change in subtropical jet strength at max shear (ms$\mathregular{^{-1}}$)", size=20)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontsize=20)
    rows = ['Winter response', 'Summer response']
    pad = 10  # in points
    for ax, row in zip(axs[:], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=24)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    plt.show()


def zonalitcz(data1, p, sigma, lat, control, csigma, clat, colormax=5):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
    data = np.copy(data1)
    fig, axs = plt.subplots(1, 1, sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Hadley Cell SNR Sensitivity", size=24, y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-colormax, colormax, 17)

    control = control / 1e11
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot = axs.contourf(meshlat, meshsigma, data, ctrs,
                        extend='both',
                        cmap=newcmap, vmin=-colormax, vmax=colormax)
    plot1 = axs.contour(cmeshlat, cmeshsigma, control, cctrs,
                        colors='k', linewidths=1.5)

    axs.invert_yaxis()
    axs.set_xlim([lat[0], lat[-1]])
    axs.set_ylim(1, 0.1)
    axs.set_xticks(np.arange(-75, 90, 15))
    axs.set_yticks(np.arange(0.2, 1.2, 0.2))
    axs.set_xticklabels(axs.get_xticks(), fontsize=20)
    axs.set_yticklabels(axs.get_yticks(), fontsize=20)
    meshsigma2 = np.copy(meshsigma)
    meshsigma1 = np.copy(meshsigma)
    # axs[i, j].contour(meshlat, meshsigma, p[i, 1-j, :], [0.005],
    # linewidths=2)
    for a in range(9):
        for b in range(34):
            if p[a, b] <= 0.005:
                meshsigma2[a, b] = np.nan
    axs.scatter(meshlat, meshsigma2, marker='o', color='gray',
                s=0.5)
    for a in range(9):
        for b in range(34):
            if p[a, b] > 0.005:
                meshsigma1[a, b] = np.nan
    axs.scatter(meshlat, meshsigma1, marker='o', color='black',
                s=5)
    plt.subplots_adjust(hspace=0.079, wspace=0.15, top=.8, bottom=.2,
                        left=.2, right=.8)
    cbar_ax = fig.add_axes([0.85, 0.2, 0.01, 0.6])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='vertical',
                     ticks=np.linspace(-colormax, colormax, 9))
    b.set_label(label='ITCZ sensitivity (deg)', size=20)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontsize=20)
    plt.show()


def zonalhad_widthsITCZ(data1, data2, p, p2, sigma, lat, control, csigma,
                        clat, colormax=4):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    import numpy as np
    data = np.copy(data1)
    fig, axs = plt.subplots(2, 2, sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Hadley Cell SNR Sensitivity", size=24, y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    cmeshlat, cmeshsigma = np.meshgrid(clat, csigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-colormax, colormax, 17)

    control = control / 1e11
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot = axs[0, 0].contourf(meshlat, meshsigma, data[1, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=-colormax, vmax=colormax)
    plot1 = axs[0, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    plot = axs[0, 1].contourf(meshlat, meshsigma, data[0, :], ctrs,
                              extend='both',
                              cmap=newcmap, vmin=-colormax, vmax=colormax)
    plot1 = axs[0, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    plot = axs[1, 0].contourf(meshlat, meshsigma, data2, ctrs,
                              extend='both',
                              cmap=newcmap, vmin=-colormax, vmax=colormax)
    plot1 = axs[1, 0].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    plot = axs[1, 1].contourf(meshlat, meshsigma, data[2, :], ctrs,
                              extend='both', cmap=newcmap, vmin=-colormax,
                              vmax=colormax)
    plot1 = axs[1, 1].contour(cmeshlat, cmeshsigma, control, cctrs,
                              colors='k', linewidths=1.5)

    axs[0, 0].set_title('ITCZ', size=20)
    axs[1, 0].set_title('Tropical belt width', size=20)
    axs[0, 1].set_title('Winter Hadley cell width', size=20)
    axs[1, 1].set_title('Summer Hadley cell width', size=20)

    meshsigma2 = np.copy(meshsigma)
    meshsigma1 = np.copy(meshsigma)
    for a in range(9):
        for b in range(34):
            if p[1, a, b] <= 0.005:
                meshsigma2[a, b] = np.nan
    axs[0, 0].scatter(meshlat, meshsigma2, marker='o', color='gray',
                      s=0.5)
    for a in range(9):
        for b in range(34):
            if p[1, a, b] > 0.005:
                meshsigma1[a, b] = np.nan
    axs[0, 0].scatter(meshlat, meshsigma1, marker='o', color='black',
                      s=5)

    meshsigma2 = np.copy(meshsigma)
    meshsigma1 = np.copy(meshsigma)
    for a in range(9):
        for b in range(34):
            if p[0, a, b] <= 0.005:
                meshsigma2[a, b] = np.nan
    axs[0, 1].scatter(meshlat, meshsigma2, marker='o', color='gray',
                      s=0.5)
    for a in range(9):
        for b in range(34):
            if p[0, a, b] > 0.005:
                meshsigma1[a, b] = np.nan
    axs[0, 1].scatter(meshlat, meshsigma1, marker='o', color='black',
                      s=5)

    meshsigma2 = np.copy(meshsigma)
    meshsigma1 = np.copy(meshsigma)
    for a in range(9):
        for b in range(34):
            if p[2, a, b] <= 0.005:
                meshsigma2[a, b] = np.nan
    axs[1, 1].scatter(meshlat, meshsigma2, marker='o', color='gray',
                      s=0.5)
    for a in range(9):
        for b in range(34):
            if p[1, a, b] > 0.005:
                meshsigma1[a, b] = np.nan
    axs[1, 1].scatter(meshlat, meshsigma1, marker='o', color='black',
                      s=5)

    meshsigma2 = np.copy(meshsigma)
    meshsigma1 = np.copy(meshsigma)
    for a in range(9):
        for b in range(34):
            if p2[a, b] <= 0.005:
                meshsigma2[a, b] = np.nan
    axs[1, 0].scatter(meshlat, meshsigma2, marker='o', color='gray',
                      s=0.5)
    for a in range(9):
        for b in range(34):
            if p2[a, b] > 0.005:
                meshsigma1[a, b] = np.nan
    axs[1, 0].scatter(meshlat, meshsigma1, marker='o', color='black',
                      s=5)

    for i in range(2):
        for j in range(2):
            axs[i, j].invert_yaxis()
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0.1)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.2, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(), fontsize=16)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(), fontsize=16)

            # axs[i, j].contour(meshlat, meshsigma, p[i, 1-j, :], [0.005],
            # linewidths=2)

    plt.subplots_adjust(hspace=0.2, wspace=0.045, top=.95, bottom=0.15,
                        left=.1, right=.9)
    cbar_ax = fig.add_axes([0.3, 0.1, 0.4, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(-colormax, colormax, 9))
    b.set_label(label='Change (deg)', size=20)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontsize=16)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    plt.show()


def zonalHad_case(strm, sigma, lat, heat, run):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    # plt.suptitle("Run " + str(run), size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    data = (strm[run, :] - strm[0, :]) / 1e10
    caxismin, caxismax, ctrs = -2, 2, np.linspace(-2, 2, 17)
    datacontrol = strm[0, :] / 1e11
    cctrs = [-1.3, -1.2, -1.1, -1, -.9, -.8, -.7, -.6, -.5, -.4, -.3, -.2, -.1,
             .1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2]

    plot2 = plt.contourf(meshlat, meshsigma, data, ctrs,
                         cmap=newcmap, vmin=caxismin, vmax=caxismax,
                         extend='both')
    m = plt.colorbar(plot2, orientation='horizontal', aspect=50,
                     format='%.1f', spacing='proportional')
    m.set_label(label="Streamfunction change (kg s$\mathregular{^{-1}}$ x 10$\mathregular{^{10}}$)\n(Contours show streamfunction and control (black and green) (kg s$\mathregular{^{-1}}$ x 10$\mathregular{^{11}}$))", fontsize=20)
    q = m.ax.get_xticklabels()
    m.ax.set_xticklabels(q, fontsize=20)
    plt.contour(meshlat, meshsigma, datacontrol, cctrs,
                colors='g', linewidths=1.5)
    plt.contour(meshlat, meshsigma, strm[run, :]/1e11, cctrs,
                colors='k', linewidths=1.5)
    ax = plt.gca()
    ax.invert_yaxis()
    if run != 0:
        sigmax = np.argmax(heat[run, :, :], axis=0)[0]
        latmax = np.argmax(heat[run, :, :], axis=1)[0]
        plt.scatter(lat[latmax], sigma[sigmax], marker='o',
                    color='deeppink', linewidth=2, s=50)
    ax.set_xlim([lat[0], lat[-1]])
    ax.set_ylim(1, 0)
    ax.set_xticks(np.arange(-75, 90, 15))
    ax.set_yticks(np.arange(0.0, 1.2, 0.2))
    ax.set_xticklabels(ax.get_xticks(), fontsize=20)
    ax.set_yticklabels(ax.get_yticks(), fontsize=20)
    plt.show()


def scatplotall(mean_diff, d_Had_vert, d_Had_itcz, gridlat):
    import matplotlib.pyplot as plt
    c1 = 'purple'
    c2 = 'red'
    c3 = 'blue'
    gridsigma = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1]
    plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    plt.axhline(gridsigma[0], linewidth=0.5, color='k')
    plt.axhline(gridsigma[1], linewidth=0.5, color='k')
    plt.axhline(gridsigma[2], linewidth=0.5, color='k')
    plt.axhline(gridsigma[3], linewidth=0.5, color='k')
    plt.axhline(gridsigma[4], linewidth=0.5, color='k')
    plt.axhline(gridsigma[5], linewidth=0.5, color='k')
    plt.axhline(gridsigma[6], linewidth=0.5, color='k')
    plt.axhline(gridsigma[7], linewidth=0.5, color='k')
    plt.axhline(gridsigma[8], linewidth=0.5, color='k')
    plt.axvline(-45.1286245, linewidth=2, color=c1)
    plt.axvline(55.88670132, linewidth=2, color=c1)
    plt.axvline(41.24746448, linewidth=2, color=c2)
    plt.axvline(18.65775303, linewidth=2, color=c3)
    plt.axvline(-30.55790252, linewidth=2, color=c2)
    plt.ylim(1.02, 0.105)
    plt.xlim(-90, 90)

    d_Had = d_Had_vert[1, :]
    d_jet = mean_diff[:, 1, :]
    d_itcz = d_Had_itcz[1]

    gridlat_2 = np.copy(gridlat)
    gridlat_0 = gridlat_2-1.2
    gridlat_1 = gridlat_2-0.6
    gridlat_3 = gridlat_2+0.6
    gridlat_4 = gridlat_2+1.2
    d_itcz = np.expand_dims(d_itcz, 0)
    deltas = np.stack((d_jet[0], d_Had[0], d_itcz[0], d_Had[1],
                       d_jet[1]), axis=0)
    gridlats = np.stack((gridlat_0, gridlat_1, gridlat_2, gridlat_3,
                         gridlat_4))
    meshy = np.zeros((9, 34*5))
    meshxx = np.ndarray.flatten(gridlats, order='F')
    meshx = np.zeros((9, 170))
    meshx[:] = meshxx

    deltas_norm = deltas*0.006
    for i in range(9):
        for j in range(34):
            meshy[i, 5*j:5*(j+1)] = gridsigma[i]-deltas_norm[:, i, j]
    plt.scatter(meshx, meshy, color=[c1, c2, c3,
                                     c2, c1], s=10)
    plt.subplots_adjust(top=.95, bottom=0.05,
                        left=.05, right=.95)


def scatplotall_trace(mean_diff, d_Had_vert, d_Had_itcz, gridlat):
    import matplotlib.pyplot as plt
    c1 = 'purple'
    c2 = 'red'
    c3 = 'blue'
    c4 = 'orange'
    c5 = 'green'
    gridsigma = [1, .9, .8, .7, .6, .5, .4, .3, .2, .1]
    plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    plt.axhline(gridsigma[0], linewidth=0.5, color='k')
    plt.axhline(gridsigma[1], linewidth=0.5, color='k')
    plt.axhline(gridsigma[2], linewidth=0.5, color='k')
    plt.axhline(gridsigma[3], linewidth=0.5, color='k')
    plt.axhline(gridsigma[4], linewidth=0.5, color='k')
    plt.axhline(gridsigma[5], linewidth=0.5, color='k')
    plt.axhline(gridsigma[6], linewidth=0.5, color='k')
    plt.axhline(gridsigma[7], linewidth=0.5, color='k')
    plt.axhline(gridsigma[8], linewidth=0.5, color='k')
    plt.axvline(-45.1286245, linewidth=2, color=c1)
    plt.axvline(55.88670132, linewidth=2, color=c2)
    plt.axvline(41.24746448, linewidth=2, color=c3)
    plt.axvline(18.65775303, linewidth=2, color=c4)
    plt.axvline(-30.55790252, linewidth=2, color=c5)
    plt.ylim(1.02, 0.105)
    plt.xlim(-90, 90)

    d_Had = d_Had_vert[1, :]
    d_jet = mean_diff[:, 1, :]
    d_itcz = d_Had_itcz[1]

    gridlat_2 = np.copy(gridlat)
    gridlat_0 = gridlat_2-1.2
    gridlat_1 = gridlat_2-0.6
    gridlat_3 = gridlat_2+0.6
    gridlat_4 = gridlat_2+1.2
    d_itcz = np.expand_dims(d_itcz, 0)
    deltas = np.stack((d_jet[0], d_Had[0], d_itcz[0], d_Had[1],
                       d_jet[1]), axis=0)
    gridlats = np.stack((gridlat_0, gridlat_1, gridlat_2, gridlat_3,
                         gridlat_4))
    meshy = np.zeros((9, 34*5))
    meshxx = np.ndarray.flatten(gridlats, order='F')
    meshx = np.zeros((9, 170))
    meshx[:] = meshxx

    deltas_norm = deltas*0.006
    for i in range(9):
        for j in range(34):
            meshy[i, 5*j:5*(j+1)] = gridsigma[i]-deltas_norm[:, i, j]
    plt.scatter(meshx, meshy, color=[c1, c5, c4,
                                     c3, c2], s=10)
    for i in range(5):
        for j in range(9):
            if i == 0:
                c6 = c1
            if i == 1:
                c6 = c5
            if i == 2:
                c6 = c4
            if i == 3:
                c6 = c3
            if i == 4:
                c6 = c2
            plt.plot(meshx[j, i::5], meshy[j, i::5], color=c6)
    plt.subplots_adjust(top=.95, bottom=0.05,
                        left=.05, right=.95)


def scatplotall_tracesub(mean_diff, d_Had_vert, d_Had_itcz, gridlat):
    import matplotlib.pyplot as plt
    c1 = 'purple'
    c2 = 'red'
    c3 = 'blue'
    c4 = 'orange'
    c5 = 'green'
    gridsigma = [.8, .6, .4]
    plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    plt.axhline(gridsigma[0], linewidth=0.5, color='k')
    plt.axhline(gridsigma[1], linewidth=0.5, color='k')
    plt.axhline(gridsigma[2], linewidth=0.5, color='k')
    plt.axvline(-45.1286245, linewidth=2, color=c1)
    plt.axvline(55.88670132, linewidth=2, color=c2)
    plt.axvline(41.24746448, linewidth=2, color=c3)
    plt.axvline(18.65775303, linewidth=2, color=c4)
    plt.axvline(-30.55790252, linewidth=2, color=c5)
    plt.ylim(0.875, 0.1)
    plt.xlim(-90, 90)

    d_Had = d_Had_vert[1, :, 2:8:2]
    d_jet = mean_diff[:, 1, 2:8:2]
    d_itcz = d_Had_itcz[1, 2:8:2]

    gridlat_2 = np.copy(gridlat)
    gridlat_0 = gridlat_2-1.2
    gridlat_1 = gridlat_2-0.6
    gridlat_3 = gridlat_2+0.6
    gridlat_4 = gridlat_2+1.2
    d_itcz = np.expand_dims(d_itcz, 0)
    deltas = np.stack((d_jet[0], d_Had[0], d_itcz[0], d_Had[1],
                       d_jet[1]), axis=0)
    gridlats = np.stack((gridlat_0, gridlat_1, gridlat_2, gridlat_3,
                         gridlat_4))
    meshy = np.zeros((3, 34*5))
    meshxx = np.ndarray.flatten(gridlats, order='F')
    meshx = np.zeros((3, 170))
    meshx[:] = meshxx

    deltas_norm = deltas*0.02
    for i in range(3):
        for j in range(34):
            meshy[i, 5*j:5*(j+1)] = gridsigma[i]-deltas_norm[:, i, j]
    plt.scatter(meshx, meshy, color=[c1, c5, c4,
                                     c3, c2], s=10)
    for i in range(5):
        for j in range(3):
            if i == 0:
                c6 = c1
            if i == 1:
                c6 = c5
            if i == 2:
                c6 = c4
            if i == 3:
                c6 = c3
            if i == 4:
                c6 = c2
            plt.plot(meshx[j, i::5], meshy[j, i::5], color=c6)
    plt.subplots_adjust(top=.95, bottom=0.05,
                        left=.05, right=.95)


def zonalhad_comp(u, duv, potemp, strm, vT, divvT, temp, z, sigma, lat, heat, run1, run2):
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    def colourscale(plotdata):
        """
        Takes data being plotted and normalises the colourscale between largest
        data value and it's negative multiple

        Parameters
        ----------
        plotdata: array
            data being plotted

        Returns
        -------
        caxismax: int
            max magnitude data value
        caxismin: int
            negative of max mag data value
        ctrs: array
            gradated colour scale contour array
        """
        M = np.nanmax(plotdata)
        m = np.nanmin(plotdata)
        if M >= abs(m):
            ctrs1 = np.arange(-M, 0, .1*M)
            ctrs2 = np.arange(0.1*M, 1.09*M, .1*M)
            ctrs = np.concatenate((ctrs1, ctrs2))
            caxismin = -M
            caxismax = M
        else:
            m = -m
            ctrs1 = np.arange(-m, 0, .1*m)
            ctrs2 = np.arange(0.1*m, 1.09*m, .1*m)
            ctrs = np.concatenate((ctrs1, ctrs2))
            caxismin = -m
            caxismax = m
        # function will not work if there exist no positive max or negative min
        return caxismin, caxismax, ctrs

    def lapserate(t, z, sigma, lat):
        """
        Produces plot of lapse rate of T data

        Parameters
        ----------
        t: array
            temperature data field
        z: array
            geopotential height of field
        sigma: array
            sigma levels of field
        lat: array
            lat of field
        """
        import numpy as np
        dT = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
        dz = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
        for i in range(np.ma.size(sigma, axis=0)-1):
            dT[i, :] = t[i+1, :] - t[i, :]
        for i in range(np.ma.size(sigma, axis=0)-1):
            dz[i, :] = z[i+1, :] - z[i, :]
        lapse = -1000 * dT[0:-1] / dz[0:-1]
        # zonalplot(lapse, sigma[0:-1], lat, 'Lapse rate')
        return lapse

    fig, axs = plt.subplots(2, 2, sharey='row', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Run " + str(run), size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    lapse = lapserate(temp[run1, :, :], z[run1, :, :], sigma, lat)
    lapse[22:, :] = 99
    lapse0 = lapserate(temp[run2, :, :], z[run2, :, :], sigma, lat)
    lapse0[22:, :] = 99

    ctrs = np.linspace(-60, 60, 17)
    ctrsa1 = np.arange(-70, 0, 7)
    ctrsb1 = np.arange(7, 70, 7)
    ctrs2 = np.concatenate((ctrsa1, ctrsb1))
    plot2 = axs[0, 0].contourf(meshlat, meshsigma, u[run1, :], ctrs,
                               cmap=newcmap, vmin=-65, vmax=65, extend='both')
    m = axs[0, 0].get_figure().colorbar(plot2, ax=axs[0, 0], pad=0.1,
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    m.set_label(label="Zonal wind (ms$\mathregular{^{-1}}$)\nContours show divergence of u'v' (1x10$\mathregular{^{-6}}$ ms$\mathregular{^{-2}}$)", fontsize=20)
    q = m.ax.get_xticklabels()
    m.ax.set_xticklabels(q, fontsize=20)
    axs[0, 0].contour(meshlat, meshsigma, duv[run2, :]*1e6, ctrs2,
                      colors='g', linewidths=1.5)
    plot1 = axs[0, 0].contour(meshlat, meshsigma, duv[run1, :]*1e6, ctrs2,
                              colors='k', linewidths=1.5)
    plt.clabel(plot1, inline=True, inline_spacing=0, fontsize=10, fmt='%.0f')

    ctrs = np.linspace(-2, 2, 17)
    #ctrs2 = [270, 280, 290, 300, 310, 320, 330, 340, 350, 400, 500]
    ctrs2 = np.linspace(200, 350, 16)
    plot3 = axs[1, 0].contourf(meshlat, meshsigma, strm[run1, :]/1e11, ctrs,
                               cmap=newcmap, vmin=-2, vmax=2, extend='both')
    n = axs[1, 0].get_figure().colorbar(plot3, ax=axs[1, 0], pad=0.1,
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    n.set_label(label="Streamfunction (1x10$\mathregular{^{11}}$ kg s$\mathregular{^{-1}}$)\nContours show potential temperature (K)", fontsize=20)
    r = n.ax.get_xticklabels()
    n.ax.set_xticklabels(r, fontsize=20)
    axs[1, 0].contour(meshlat, meshsigma, potemp[run2, :], ctrs2,
                      colors='g', linewidths=1.5)
    plot1 = axs[1, 0].contour(meshlat, meshsigma, potemp[run1, :], ctrs2,
                              colors='k', linewidths=1.5)
    axs[1, 0].clabel(plot1, inline=True, inline_spacing=-10, fontsize=10,
                     fmt='%.1f')

    cmax = 6
    ctrs = np.linspace(-cmax, cmax, 17)
    cmin = -cmax
    ctrsc1 = np.arange(-30, 0, 3)
    ctrsd1 = np.arange(3, 30, 3)
    ctrs2 = np.concatenate((ctrsc1, ctrsd1))
    plot4 = axs[0, 1].contourf(meshlat, meshsigma, u[run1, :]-u[run2, :], ctrs,
                               cmap=newcmap, vmin=cmin, vmax=cmax,
                               extend='both')
    o = axs[0, 1].get_figure().colorbar(plot4, ax=axs[0, 1], pad=0.1,
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    o.set_label(label="Zonal wind change (ms$\mathregular{^{-1}}$)\nContours show v'T' (K ms$\mathregular{^{-1}}$)",
                fontsize=20)
    s = o.ax.get_xticklabels()
    o.ax.set_xticklabels(s, fontsize=20)
    axs[0, 1].contour(meshlat, meshsigma, vT[run2, :], ctrs2, colors='g',
                      linewidths=1.5)
    plot1 = axs[0, 1].contour(meshlat, meshsigma, vT[run1, :], ctrs2,
                              colors='k', linewidths=1.5)
    axs[0, 1].clabel(plot1, inline=True, inline_spacing=0, fontsize=10,
                     fmt='%.0f')

    cmax = 2
    ctrs = np.linspace(-cmax, cmax, 17)
    cmin = -cmax
    ctrsc = np.arange(-21, 0, 3)
    ctrsd = np.arange(3, 24, 3)
    ctrs2 = np.concatenate((ctrsc, ctrsd))
    plot5 = axs[1, 1].contourf(meshlat, meshsigma,
                               (strm[run1, :]-strm[run2, :])/1e10, ctrs,
                               cmap=newcmap, vmin=cmin, vmax=cmax,
                               extend='both')
    p = axs[1, 1].get_figure().colorbar(plot5, ax=axs[1, 1], pad=0.1,
                                        orientation='horizontal', aspect=50,
                                        format='%.1f', spacing='proportional')
    t = p.ax.get_xticklabels()
    p.ax.set_xticklabels(t, fontsize=20)
    p.set_label(label="Streamfunction change (1x10$\mathregular{^{10}}$ kg s$\mathregular{^{-1}}$)\nContours show divergence of v'T' (1x10$\mathregular{^{-6}}$ K s$\mathregular{^{-1}}$)", fontsize=20)    
    axs[1, 1].contour(meshlat, meshsigma, divvT[run2, :]*1e6, ctrs2,
                      colors='g', linewidths=1.5)
    plot1 = axs[1, 1].contour(meshlat, meshsigma, divvT[run1, :]*1e6, ctrs2,
                              colors='k', linewidths=1.5)
    axs[1, 1].clabel(plot1, inline=True, inline_spacing=0, fontsize=10,
                     fmt='%.0f')
    for i in range(2):
        for j in range(2):
            axs[i, j].contour(meshlat[1:, :], meshsigma[1:, :], lapse0, [2], colors='darkblue')
            axs[i, j].contour(meshlat[1:, :], meshsigma[1:, :], lapse, [2],
                              linewidths=2, colors='darkblue')
            axs[i, j].invert_yaxis()
            if run1 != 0:
                sigmax = np.argmax(heat[run1, :, :], axis=0)[0]
                latmax = np.argmax(heat[run1, :, :], axis=1)[0]
                axs[i, j].scatter(lat[latmax], sigma[sigmax], marker='o',
                                  color='deeppink', linewidth=2, s=50)
            if run2 != 0:
                sigmax = np.argmax(heat[run2, :, :], axis=0)[0]
                latmax = np.argmax(heat[run2, :, :], axis=1)[0]
                axs[i, j].scatter(lat[latmax], sigma[sigmax], marker='o',
                                  color='g', linewidth=2, s=50)
            axs[i, j].set_xlim([lat[0], lat[-1]])
            axs[i, j].set_ylim(1, 0)
            axs[i, j].set_xticks(np.arange(-75, 90, 15))
            axs[i, j].set_yticks(np.arange(0.0, 1.2, 0.2))
            axs[i, j].set_xticklabels(axs[i, j].get_xticks(),
                                      fontsize=20)
            axs[i, j].set_yticklabels(axs[i, j].get_yticks(),
                                      fontsize=20)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=24, fontweight='bold')
    plt.subplots_adjust(hspace=0.2, wspace=.05, top=.97, bottom=0.1, left=.05,
                        right=.95)
    plt.show()


def zonalhad_temp(potemp, temp, sigma, lat, heat, run1, run2):
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    fig, axs = plt.subplots(1, 1, facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Run " + str(run), size=30, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)

    ctrs = np.linspace(-.02, .02, 17)
    ctrs2 = np.linspace(270, 350, 9)
    plot3 = axs.contourf(meshlat, meshsigma, np.log(temp[run1]/temp[run2]), ctrs,
                         cmap=newcmap, extend='both')
    n = axs.get_figure().colorbar(plot3, ax=axs, pad=0.1,
                                  orientation='horizontal', aspect=50,
                                  format='%.1f', spacing='proportional')
    n.set_label(label="Temperature difference (K)\nContours show potential temperature (K)", fontsize=20)
    r = n.ax.get_xticklabels()
    n.ax.set_xticklabels(r, fontsize=20)
    axs.contour(meshlat, meshsigma, potemp[run2, :], ctrs2,
                colors='g', linewidths=1.5)
    plot1 = axs.contour(meshlat, meshsigma, potemp[run1, :], ctrs2,
                        colors='k', linewidths=1.5)
    axs.clabel(plot1, inline=True, inline_spacing=-10, fontsize=10, fmt='%.1f')

    axs.invert_yaxis()
    if run1 != 0:
        sigmax = np.argmax(heat[run1, :, :], axis=0)[0]
        latmax = np.argmax(heat[run1, :, :], axis=1)[0]
        axs.scatter(lat[latmax], sigma[sigmax], marker='o',
                    color='deeppink', linewidth=2, s=50)
    if run2 != 0:
        sigmax = np.argmax(heat[run2, :, :], axis=0)[0]
        latmax = np.argmax(heat[run2, :, :], axis=1)[0]
        axs.scatter(lat[latmax], sigma[sigmax], marker='o',
                    color='g', linewidth=2, s=50)
    axs.set_xlim([lat[0], lat[-1]])
    axs.set_ylim(1, 0)
    axs.set_xticks(np.arange(-75, 90, 15))
    axs.set_yticks(np.arange(0.0, 1.2, 0.2))
    axs.set_xticklabels(axs.get_xticks(),
                        fontsize=20)
    axs.set_yticklabels(axs.get_yticks(),
                        fontsize=20)
    plt.subplots_adjust(hspace=0.2, wspace=.05, top=.97, bottom=0.1, left=.05,
                        right=.95)
    plt.show()


def zonalhad_rean(data_id, data_djf, data_jja, u_id, u_djf, u_jja, sigma,
                  lat, p, lat_era):
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    data_id = np.copy(data_id)
    data_djf = np.copy(data_djf)
    data_jja = np.copy(data_jja)
    fig, axs = plt.subplots(3, 1, sharex='col', facecolor='w',
                            edgecolor='k', linewidth=2)
    # plt.suptitle("Hadley Cell Sensitivity", size=24, y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    meshlat_era, meshp = np.meshgrid(lat_era, p)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[238:274, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    data_id = data_id / 1e11
    data_djf = data_djf / 1e11
    data_jja = data_jja / 1e11
    caxismin, caxismax, ctrs = -2, 2, np.linspace(-2, 2, 21)
    caxismin1, caxismax1, ctr1 = -50, 50, np.array([-50,-45,-40,-35,-30,-25,-20,-15,-10,-5,5,10,15,20,25,30,35,40,45,50])

    plot = axs[0].contourf(meshlat, meshsigma, data_id, ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    axs[0].contour(meshlat, meshsigma, u_id, ctr1,
                      colors='k', linewidths=1.5)
    #  axs[0, 0].clabel(plot1, [-1,-.5,.5,1], inline=True, inline_spacing=-3,
    #                   fontsize=10,
    #                   fmt='%.1f')

    plot = axs[1].contourf(meshlat_era, meshp, data_djf, ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    axs[1].contour(meshlat_era, meshp, u_djf, ctr1,
                      colors='k', linewidths=1.5)
    plot = axs[2].contourf(meshlat_era, meshp, data_jja, ctrs,
                              extend='both',
                              cmap=newcmap, vmin=caxismin, vmax=caxismax)
    axs[2].contour(meshlat_era, meshp, u_jja, ctr1,
                      colors='k', linewidths=1.5)

    for i in range(3):
        axs[i].invert_yaxis()
        axs[i].set_xlim([lat[0], lat[-1]])
        axs[i].set_xticks(np.arange(-75, 90, 15))
        axs[i].set_xticklabels(axs[i].get_xticks(), fontsize=16)

    axs[0].set_ylim(1, 0)
    axs[0].set_yticks(np.arange(0.2, 1.2, 0.2))
    axs[0].set_yticklabels(axs[0].get_yticks(), fontsize=16)
    axs[1].set_ylim(1000, 0)
    axs[1].set_yticks(np.arange(200, 1200, 200))
    axs[1].set_yticklabels(axs[1].get_yticks(), fontsize=16)
    axs[2].set_ylim(1000, 0)
    axs[2].set_yticks(np.arange(200, 1200, 200))
    axs[2].set_yticklabels(axs[2].get_yticks(), fontsize=16)

    plt.subplots_adjust(hspace=0.08, wspace=0.045, top=.95, bottom=0.15,
                        left=.3, right=.7)
    cbar_ax = fig.add_axes([0.3, 0.1, 0.4, 0.015])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='uniform', extend='max',
                     orientation='horizontal',
                     ticks=np.linspace(caxismin, caxismax, 9))
    b.set_label(label='Streamfunction (kg s$\mathregular{^{-1}}$ x 10$\mathregular{^{11}}$)', size=20)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontsize=16)

    rows = ['Idealized GCM', 'DJF', 'JJA']
    pad = 10  # in points
    for ax, row in zip(axs[:], rows):
        ax.annotate(row, xy=(0, 0), xytext=(-ax.yaxis.labelpad - pad, 0),
                    xycoords=ax.yaxis.label, textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=24)
    for ax, lab in zip(axs, ['a', 'b', 'c']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=18, fontweight='bold')
    plt.show()
