#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  2 11:03:09 2018

@author: bakerh
"""
import numpy as np


def mapplot60(plotdata_sig, plotdata, mask='yes', title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from netcdfread import ncread
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 1.875)
    long = np.concatenate((lon, [360]))
    fig, axs = plt.subplots(2, 2, facecolor='w',
                            edgecolor='k', linewidth=2)
    # fig.suptitle('North Atlantic Jet Sensitivity', fontsize=20, y=0.96)
    for a in range(2):
        for b in range(2):
            if b == 0:
                mx = 1.2
            else:
                mx = .8
            plotdata1 = np.concatenate((plotdata_sig[a, 1-b, :],
                                        np.expand_dims(plotdata_sig[a, 1-b, :, 0],
                                        axis=1)), axis=1) * 1e5
            plotdata3 = np.concatenate((plotdata[a, 1-b, :],
                                        np.expand_dims(plotdata[a, 1-b, :, 0],
                                        axis=1)), axis=1) * 1e5
            if mask == 'yes':
                lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
                lsm = np.concatenate((lsm, np.expand_dims(lsm[:, 0], axis=1)),
                                     axis=1)
                plotdata1 = np.ma.masked_array(plotdata1, lsm)
                plotdata3 = np.ma.masked_array(plotdata3, lsm)
            meshlon, meshlat = np.meshgrid(long, lat)
            ax1 = axs[a, b]
            m = Basemap(projection='cyl', llcrnrlat=-60, urcrnrlat=60,
                        llcrnrlon=0, urcrnrlon=360, resolution='c', ax=ax1)
            m.drawcoastlines()
            m.drawmapboundary(linewidth=2)
            x, y = m(meshlon, meshlat)
            mycmap2 = plt.cm.YlOrRd(np.arange(256))
            mycmap1 = plt.cm.Blues_r(np.arange(256))
            my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
            #my_cmap[239:274, :] = 1
            newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet",
                                                                   my_cmap)
            caxismin, caxismax, ctrs = -mx, mx, np.linspace(-mx, mx, 17)
            plot = m.contourf(x, y, plotdata1, ctrs,
                              cmap=newcmap, vmin=caxismin, vmax=caxismax,
                              extend='both')
            m.contour(x, y, plotdata3, ctrs, colors='k')
            ax1.set_ylim(-60, 60)
            m.drawparallels(np.arange(-60, 90, 30),
                            labels=[True, False, False, True], linewidth=0)
            m.drawmeridians(np.arange(0, 390, 30),
                            labels=[True, False, False, True], linewidth=0)
            ax1.set_yticks(np.arange(-60., 90., 30.))
            ax1.set_xticks(np.arange(0., 390., 30.))
            ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                            direction='out', length=5, width=2)
            ax1.set_yticks(np.arange(-60., 70., 10.), minor=True)
            ax1.set_xticks(np.arange(0., 370., 10.), minor=True)
            ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                            direction='out', length=4, width=1)
            if a == 0 and b == 0:
                ax1.set_ylabel('DJF', fontsize=16, labelpad=25)
                ax1.set_title('Latitude', fontsize=16, y=1.08)
            if a == 1 and b == 0:
                c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
                ax1.set_ylabel('JJA', fontsize=16, labelpad=25)
                c.set_label(label='Poleward jet latitude shift ($^\circ$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
            if a == 0 and b == 1:
                ax1.set_title('Speed', fontsize=16, y=1.08)
            if a == 1 and b == 1:
                c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
                c.set_label(label='Jet speed increase (ms$\mathregular{^{-1}}$ [$\mathregular{1x10^5km^2}$K]$\mathregular{^{-1}}$)', size=16)
    plt.subplots_adjust(hspace=0, wspace=0.1, top=.95, bottom=0.05, left=.05,
                        right=.95)

    plt.show()
