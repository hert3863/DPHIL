#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 12:00:32 2018

@author: bakerh
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  5 13:45:17 2018

@author: bakerh
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netcdfread import ncread
from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic


def ept():
    from netcdfread import ncread
    import glob
    from netCDF4 import Dataset

    def gradient(u, lat, lon):
        meshlon, meshlat = np.meshgrid(lon, lat)
        phi = meshlat*np.pi/180
        cs_phi = np.cos(phi)
        dphi = lat*np.pi/180
        dtheta = lon*np.pi/180
        a = 6.371e6

        u_y = np.gradient(u, dphi, axis=1) / a
        u_x = np.gradient(u, dtheta, axis=2) / (a*cs_phi)
        return u_x, u_y

    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    hgt = ncread('orog.nc', 'ht')[0, 0]

    sindices = np.zeros((10))
    for i in range(10):
        sindices[i] = i*12
    sindices = sindices.astype(int)

    ThetaE = {}
    ThetaE_x = {}
    ThetaE_y = {}
    exps = ['Plus15-Future_HCO2', 'Plus15-Future_LCO2']
    exps = ['All-Nat', 'SST-Nat', 'GHG-Nat']

    for exp in exps:
        a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                      exp + '/mon/hurs/*')
        t_e_all = np.zeros((12, 145, 192))
        t_ex_all = np.zeros((12, 145, 192))
        t_ey_all = np.zeros((12, 145, 192))
        for g, d in enumerate(a):
            b = d[:60]+str('psl')+d[64:69]+str(16222)+d[73:]
            z = d[:60]+str('tas')+d[64:69]+str(3236)+d[73:]

            nc_fid = Dataset(d, 'r')
            nc_fib = Dataset(b, 'r')
            nc_fiz = Dataset(z, 'r')

            T = nc_fiz.variables['item3236_monthly_mean'
                                 ][:, 0, :]
            rh = nc_fid.variables['item3245_monthly_mean'
                                  ][:, 0, :]
            mslp = nc_fib.variables['item16222_monthly_mean'
                                    ][:, 0, :]

            p1 = mslp*np.exp(-hgt/(29.27*T))/100
            E = 6.1094*np.exp(17.625*(T-273.15)/(T-30.11))
            e = rh*E/100
            r = 0.622 * e / (p1-e)
            Td = 273.15 + 243.04*(np.log(rh/100)+((17.625*(T-273.15))/(T-30.11)))/(17.625-np.log(rh/100)-((17.625*(T-273.15))/(T-30.11)))
            Tl = 56 + 1/(1/(Td-56)+np.log(T/Td)/800)
            ThetaL = T*((1000/(p1-e))**0.2854) * (T/Tl)**(0.28*r)
            t_e = ThetaL * np.exp((3036/Tl-1.78)*r*(1+0.448*r))
            '''
            p0 = 1000
            p1 = mslp*np.exp(-hgt/(29.27*T))/100
            Lv = 2.51e6  # latent heat of vap at the triple point [J/kg] Bluestein,p203
            cpd = 1004.64  # specific heat at constant pressure for air
            R = 287.04  # specific gas constant for air [J/(Kg-K)]
            kap = R/cpd

            # r = mixhum_ptrh( p*0.01, t, h, 1 ) # convert relative humidity to mixing ratio
            E = 6.1094*np.exp(17.625*(T-273.15)/(T-30.11))
            e = rh*E/100
            r = 0.622 * e / (p1-e)
            t_e = T+(Lv/cpd)*r
            t_e = t_e*(p0/p1)**kap
            '''
            t_ex, t_ey = gradient(t_e, lat, lon)
            for i in range(12):
                t_e_all[i] += np.mean(t_e[sindices+i], axis=(0))
                t_ex_all[i] += np.mean(t_ex[sindices+i], axis=(0))
                t_ey_all[i] += np.mean(t_ey[sindices+i], axis=(0))
            print('Done ' + exp + ' ' + str(g+1))
        ThetaE[exp] = t_e_all/len(a)
        ThetaE_x[exp] = t_ex_all/len(a)
        ThetaE_y[exp] = t_ey_all/len(a)

    return ThetaE, ThetaE_x, ThetaE_y


def potemp1(t, p, h):

    '''
    ; Calculate Equivalent Potential Temperature
    ; ...  https://en.wikipedia.org/wiki/Equivalent_potential_temperature
    ;
    ;  p : pressure (must be scalar, 1D or same size as t & q)
    ;  t : temperature (degK)
    ;  h : humidity variable
    ;      (a) specific humidity (kg/kg) : humVarType = "h"
    ;      (b) mixing ratio (kg/kg/)     : humVarType = "r"
    ;      (c) dew point (K)             : humVarType = "td"
    ;      (d) relative humidity (%)     : humVarType = "rh"
    ;
    ;  p0: reference pressure: must be same units as p (Pa)
    '''
    p0 = 100000
    Lv = 2.51e6  # latent heat of vap at the triple point [J/kg] Bluestein,p203
    cpd = 1004.64  # specific heat at constant pressure for air
    R = 287.04  # specific gas constant for air [J/(Kg-K)]
    kap = R/cpd

    # r = mixhum_ptrh( p*0.01, t, h, 1 ) # convert relative humidity to mixing ratio
    E = 6.1094*np.exp(17.625*(t-273.15)/(t-30.11))
    e = h*E/100
    r = 0.622 * e / (0.01*p-e)
    t_e = t+(Lv/cpd)*r
    t_e = t_e*(p0/p)**kap
    return t_e


def fig1(had_l, had_h, mir_l, mir_h, cam_l, cam_h,
         had_lp, had_hp, mir_lp, mir_hp, cam_lp, cam_hp):
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lat_cam = np.linspace(-90, 90, 96)
    lon_cam = np.linspace(0, 357.5, 144)
    lat_mir = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon_mir = np.arange(0, 360, 360/256)
    h_data = np.mean(had_h[5:8], axis=0) - np.mean(had_l[5:8], axis=0)
    m_data = np.mean(mir_h[5:8], axis=0) - np.mean(mir_l[5:8], axis=0)
    c_data = (np.mean(cam_h[5:8], axis=0) - np.mean(cam_l[5:8], axis=0))*2.5
    h_datap = np.mean(had_hp[5:8], axis=0) - np.mean(had_lp[5:8], axis=0)
    m_datap = np.mean(mir_hp[5:8], axis=0) - np.mean(mir_lp[5:8], axis=0)
    c_datap = (np.mean(cam_hp[5:8], axis=0) - np.mean(cam_lp[5:8], axis=0))*2.5

    h_data, lon1 = shiftgrid(180., h_data, lon, start=False)
    h_data, lon1 = addcyclic(h_data, lon1)
    m_data, lon1_mir = shiftgrid(180., m_data, lon_mir, start=False)
    m_data, lon1_mir = addcyclic(m_data, lon1_mir)
    c_data, lon1_cam = shiftgrid(180., c_data, lon_cam, start=False)
    c_data, lon1_cam = addcyclic(c_data, lon1_cam)
    h_datap, lon1 = shiftgrid(180., h_datap, lon, start=False)
    h_datap, lon1 = addcyclic(h_datap, lon1)
    m_datap, lon1_mir = shiftgrid(180., m_datap, lon_mir, start=False)
    m_datap, lon1_mir = addcyclic(m_datap, lon1_mir)
    c_datap, lon1_cam = shiftgrid(180., c_datap, lon_cam, start=False)
    c_datap, lon1_cam = addcyclic(c_datap, lon1_cam)

    meshlon, meshlat = np.meshgrid(lon1, lat)
    meshlon_mir, meshlat_mir = np.meshgrid(lon1_mir, lat_mir)
    meshlon_cam, meshlat_cam = np.meshgrid(lon1_cam, lat_cam)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-1, 1, 17)
    #ctrs2 = np.array([-1, -.875, -.75, -.625, -.5, -.375, -.25, -.125, .125, .25, .375, .5, .625, .75, .875, 1])
    ctrs2 = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1])
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(3, 1, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, h_datap, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contour(x, y, h_data, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax2 = fig.add_subplot(3, 1, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon_mir, meshlat_mir)
    m.contourf(x, y, m_datap, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contour(x, y, m_data, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax3 = fig.add_subplot(3, 1, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon_cam, meshlat_cam)
    plot = m.contourf(x, y, c_datap, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                      vmax=np.max(ctrs), extend='both')
    m.contour(x, y, c_data, ctrs2, colors='k')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    cbar_ax = fig.add_axes([0.4, 0.04, 0.2, 0.005])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=12,
                fontsize=12, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 40  # in points

    ax1.annotate('HadAM3P', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax2.annotate('MIROC5', xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax2.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('CAM4', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)

    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=0, top=.99, bottom=0.085,
                        left=0, right=1)


def fig2(nat, sst, ghg, amipco2, amip4K, natp, sstp, ghgp, amipco2p, amip4Kp,
         te_nat, te_sst, te_ghg):
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')

    g_data = np.mean(ghg[5:8], axis=0) - np.mean(nat[5:8], axis=0)
    s_data = np.mean(sst[5:8], axis=0) - np.mean(nat[5:8], axis=0)
    g_datap = np.mean(ghgp[5:8], axis=0) - np.mean(natp[5:8], axis=0)
    s_datap = np.mean(sstp[5:8], axis=0) - np.mean(natp[5:8], axis=0)
    teg_data = np.mean(te_ghg[5:8], axis=0) - np.mean(te_nat[5:8], axis=0)
    tes_data = np.mean(te_sst[5:8], axis=0) - np.mean(te_nat[5:8], axis=0)

    g_data, lon1 = shiftgrid(180., g_data, lon, start=False)
    g_data, lon1 = addcyclic(g_data, lon1)
    s_data, lon1 = shiftgrid(180., s_data, lon, start=False)
    s_data, lon1 = addcyclic(s_data, lon1)
    g_datap, lon1 = shiftgrid(180., g_datap, lon, start=False)
    g_datap, lon1 = addcyclic(g_datap, lon1)
    s_datap, lon1 = shiftgrid(180., s_datap, lon, start=False)
    s_datap, lon1 = addcyclic(s_datap, lon1)
    teg_data, lon1 = shiftgrid(180., teg_data, lon, start=False)
    teg_data, lon1 = addcyclic(teg_data, lon1)
    tes_data, lon1 = shiftgrid(180., tes_data, lon, start=False)
    tes_data, lon1 = addcyclic(tes_data, lon1)
    a4_data, lon1 = shiftgrid(180., amip4K.mean(axis=0), lon, start=False)
    a4_data, lon1 = addcyclic(a4_data, lon1)
    ac_data, lon1 = shiftgrid(180., amipco2.mean(axis=0), lon, start=False)
    ac_data, lon1 = addcyclic(ac_data, lon1)
    a4p_data, lon1 = shiftgrid(180., amip4Kp.mean(axis=0), lon, start=False)
    a4p_data, lon1 = addcyclic(a4p_data, lon1)
    acp_data, lon1 = shiftgrid(180., amipco2p.mean(axis=0), lon, start=False)
    acp_data, lon1 = addcyclic(acp_data, lon1)

    meshlon, meshlat = np.meshgrid(lon1, lat)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-1, 1, 17)
    #ctrs2 = np.array([-1, -.875, -.75, -.625, -.5, -.375, -.25, -.125, .125, .25, .375, .5, .625, .75, .875, 1])
    ctrs2 = np.array([-2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6])
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap1 = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap1)
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(3, 2, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, acp_data, ctrs*3, cmap=newcmap, vmin=np.min(ctrs*3),
               vmax=np.max(ctrs*3), extend='both')
    m.contour(x, y, ac_data, ctrs2*2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax2 = fig.add_subplot(3, 2, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot1 = m.contourf(x, y, a4p_data, ctrs*3, cmap=newcmap, vmin=np.min(ctrs*3),
                       vmax=np.max(ctrs*3), extend='both')
    m.contour(x, y, a4_data, ctrs2*2, colors='k')
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax3 = fig.add_subplot(3, 2, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, g_datap, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contour(x, y, g_data, ctrs2, colors='k')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax4 = fig.add_subplot(3, 2, 4)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot2 = m.contourf(x, y, s_datap, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                       vmax=np.max(ctrs), extend='both')
    m.contour(x, y, s_data, ctrs2, colors='k')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    ax4.set_yticks(np.arange(-30., 120., 30.))
    ax4.set_xticks(np.arange(-180., 270., 90.))
    ax4.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax4.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax4.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax4.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    cbar_ax = fig.add_axes([0.94, 0.70, 0.0075, 0.25])
    b = fig.colorbar(plot1, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16) 
    cbar_ax = fig.add_axes([0.94, 0.375, 0.0075, 0.25])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16) 

    ax1.annotate('AMIP4xCO2', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('AMIP4K', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax3.annotate('GHG', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax4.annotate('SST', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=0.05, top=.97, bottom=0.03,
                        left=0.04, right=0.925)


def fig2b(nat, sst, ghg, amipco2, amip4K, natp, sstp, ghgp, amipco2p, amip4Kp,
         te_nat, te_sst, te_ghg):
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')

    g_data = np.mean(ghg[5:8], axis=0) - np.mean(nat[5:8], axis=0)
    s_data = np.mean(sst[5:8], axis=0) - np.mean(nat[5:8], axis=0)
    g_datap = np.mean(ghgp[5:8], axis=0) - np.mean(natp[5:8], axis=0)
    s_datap = np.mean(sstp[5:8], axis=0) - np.mean(natp[5:8], axis=0)
    teg_data = np.mean(te_ghg[5:8], axis=0) - np.mean(te_nat[5:8], axis=0)
    tes_data = np.mean(te_sst[5:8], axis=0) - np.mean(te_nat[5:8], axis=0)

    g_data, lon1 = shiftgrid(180., g_data, lon, start=False)
    g_data, lon1 = addcyclic(g_data, lon1)
    s_data, lon1 = shiftgrid(180., s_data, lon, start=False)
    s_data, lon1 = addcyclic(s_data, lon1)
    g_datap, lon1 = shiftgrid(180., g_datap, lon, start=False)
    g_datap, lon1 = addcyclic(g_datap, lon1)
    s_datap, lon1 = shiftgrid(180., s_datap, lon, start=False)
    s_datap, lon1 = addcyclic(s_datap, lon1)
    teg_data, lon1 = shiftgrid(180., teg_data, lon, start=False)
    teg_data, lon1 = addcyclic(teg_data, lon1)
    tes_data, lon1 = shiftgrid(180., tes_data, lon, start=False)
    tes_data, lon1 = addcyclic(tes_data, lon1)
    a4_data, lon1 = shiftgrid(180., amip4K, lon, start=False)
    a4_data, lon1 = addcyclic(a4_data, lon1)
    ac_data, lon1 = shiftgrid(180., amipco2, lon, start=False)
    ac_data, lon1 = addcyclic(ac_data, lon1)
    a4p_data, lon1 = shiftgrid(180., amip4Kp, lon, start=False)
    a4p_data, lon1 = addcyclic(a4p_data, lon1)
    acp_data, lon1 = shiftgrid(180., amipco2p, lon, start=False)
    acp_data, lon1 = addcyclic(acp_data, lon1)

    meshlon, meshlat = np.meshgrid(lon1, lat)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-1, 1, 17)*.8
    #ctrs2 = np.array([-1, -.875, -.75, -.625, -.5, -.375, -.25, -.125, .125, .25, .375, .5, .625, .75, .875, 1])
    #ctrs = np.array([-5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -.5, -.4, -.3, -.2, -.1, .1, .2, .3, .4, .5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0])/5
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    #my_cmap[230:282, :] = 1
    my_cmap1 = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap1)
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax5 = fig.add_subplot(3, 2, 5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, teg_data, ctrs*5, cmap=newcmap1, vmin=np.min(ctrs*5),
               vmax=np.max(ctrs*5), extend='both')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax5.set_yticks(np.arange(-30., 120., 30.))
    ax5.set_xticks(np.arange(-180., 270., 90.))
    ax5.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax5.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax5.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax5.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax5.set_title('a', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax6 = fig.add_subplot(3, 2, 6)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot3 = m.contourf(x, y, tes_data, ctrs*5, cmap=newcmap1, vmin=np.min(ctrs*5),
                       vmax=np.max(ctrs*5), extend='both')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    ax6.set_yticks(np.arange(-30., 120., 30.))
    ax6.set_xticks(np.arange(-180., 270., 90.))
    ax6.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax6.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax6.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax6.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax6.set_title('b', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    cbar_ax = fig.add_axes([0.94, 0.051, 0.0075, 0.25])
    b = fig.colorbar(plot3, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (K)', size=16,
                fontsize=16, fontname='Arial', labelpad=1)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16) 

    ax5.annotate(r'GHG $\theta_e$', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax6.annotate(r'SST $\theta_e$', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    plt.show()
    plt.subplots_adjust(hspace=0.12, wspace=0.05, top=.97, bottom=0.03,
                        left=0.04, right=0.925)


def fig3(reg_had, reg_mir, reg_cam, reg_ghg, reg_sst):
    '''
    mpl.rcParams['hatch.linewidth'] = 1.0
    def master_reg(v_l, v_h, t_l, t_h, mod='had'):
        from scipy import stats
        v_anom = v_h - v_l.mean(axis=0)
        var_anom = t_h - t_l.mean(axis=0)

        if mod == 'had':
            la = 145
            ll = 192
            circ_ind = (v_anom[:, 30, 0] - v_anom[:, 32, 21] +
                        v_anom[:, 31, 42] - v_anom[:, 37, 94] +
                        v_anom[:, 39, 121] - v_anom[:, 34, 138] +
                        v_anom[:, 35, 154] - v_anom[:, 30, 172])
            lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        elif mod == 'mir':
            la = 128
            ll = 256
            circ_ind = (v_anom[:, 95, 252] - v_anom[:, 95, 34] +
                        v_anom[:, 106, 73] - v_anom[:, 97, 145] +
                        v_anom[:, 97, 165] - v_anom[:, 98, 185] +
                        v_anom[:, 99, 205] - v_anom[:, 102, 232])
            lsm = ncread('/home/bakerh/Documents/DPhil/Python/lsm_t85.nc', 'lsm')
        elif mod == 'cam':
            la = 96
            ll = 144
            circ_ind = (v_anom[:, 73, 143] - v_anom[:, 72, 17] +
                        v_anom[:, 72, 32] + v_anom[:, 71, 92] -
                        v_anom[:, 71, 105] +
                        v_anom[:, 71, 116] - v_anom[:, 71, 127])
            lsm = ncread('/home/bakerh/Documents/DPhil/Python/lsm_n72.nc', 'lsm')
        circ_ind = (circ_ind - np.mean(circ_ind)) / circ_ind.std()
        mask = np.ones(np.shape(lsm))
        mask *= lsm
        mask /= np.max(mask)
        beta = np.zeros((la, ll))
        sig = np.zeros((la, ll))
        for i in range(la):
            for j in range(ll):
                hold = stats.linregress(circ_ind[:], (var_anom[:, i, j]))
                beta[i, j] = hold[0]
                sig[i, j] = ((hold[3])<=0.05).astype(int) * mask[i, j]
        return beta, sig
    '''

    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lat_cam = np.linspace(-90, 90, 96)
    lon_cam = np.linspace(0, 357.5, 144)
    lat_mir = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon_mir = np.arange(0, 360, 360/256)

    h_data = np.copy(reg_had)[0]
    m_data = np.copy(reg_mir)[0]
    c_data = np.copy(reg_cam)[0]
    g_data = np.copy(reg_ghg)[0]
    s_data = np.copy(reg_sst)[0]
    h_datas = np.copy(reg_had)[1]
    m_datas = np.copy(reg_mir)[1]
    c_datas = np.copy(reg_cam)[1]
    g_datas = np.copy(reg_ghg)[1]
    s_datas = np.copy(reg_sst)[1]

    h_data, lon1 = shiftgrid(180., h_data, lon, start=False)
    h_data, lon1 = addcyclic(h_data, lon1)
    m_data, lon1_mir = shiftgrid(180., m_data, lon_mir, start=False)
    m_data, lon1_mir = addcyclic(m_data, lon1_mir)
    c_data, lon1_cam = shiftgrid(180., c_data, lon_cam, start=False)
    c_data, lon1_cam = addcyclic(c_data, lon1_cam)
    g_data, lon1 = shiftgrid(180., g_data, lon, start=False)
    g_data, lon1 = addcyclic(g_data, lon1)
    s_data, lon1 = shiftgrid(180., s_data, lon, start=False)
    s_data, lon1 = addcyclic(s_data, lon1)
    h_datas, lon1 = shiftgrid(180., h_datas, lon, start=False)
    h_datas, lon1 = addcyclic(h_datas, lon1)
    m_datas, lon1_mir = shiftgrid(180., m_datas, lon_mir, start=False)
    m_datas, lon1_mir = addcyclic(m_datas, lon1_mir)
    c_datas, lon1_cam = shiftgrid(180., c_datas, lon_cam, start=False)
    c_datas, lon1_cam = addcyclic(c_datas, lon1_cam)
    g_datas, lon1 = shiftgrid(180., g_datas, lon, start=False)
    g_datas, lon1 = addcyclic(g_datas, lon1)
    s_datas, lon1 = shiftgrid(180., s_datas, lon, start=False)
    s_datas, lon1 = addcyclic(s_datas, lon1)

    meshlon, meshlat = np.meshgrid(lon1, lat)
    meshlon_mir, meshlat_mir = np.meshgrid(lon1_mir, lat_mir)
    meshlon_cam, meshlat_cam = np.meshgrid(lon1_cam, lat_cam)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-.6, .6, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(5, 1, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, h_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contourf(x, y, h_datas, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.85)
    m.scatter(x[30, 0+96], y[30, 0+96], color='fuchsia', s=25)
    m.scatter(x[32, 21+96], y[32, 21+96], color='green', s=25)
    m.scatter(x[31, 42+96], y[31, 42+96], color='fuchsia', s=25)
    m.scatter(x[37, 94-96], y[37, 94-96], color='green', s=25)
    m.scatter(x[39, 121-96], y[39, 121-96], color='fuchsia', s=25)
    m.scatter(x[34, 138-96], y[34, 138-96], color='green', s=25)
    m.scatter(x[35, 154-96], y[35, 154-96], color='fuchsia', s=25)
    m.scatter(x[30, 172-96], y[30, 172-96], color='green', s=25)

    ax2 = fig.add_subplot(5, 1, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon_mir, meshlat_mir)
    m.contourf(x, y, m_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contourf(x, y, m_datas, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.85)
    m.scatter(x[95, 252-128], y[95, 252-128], color='fuchsia', s=25)
    m.scatter(x[95, 34+128], y[95, 34+128], color='green', s=25)
    m.scatter(x[106, 73+128], y[106, 73+128], color='fuchsia', s=25)
    m.scatter(x[97, 145-128], y[97, 145-128], color='green', s=25)
    m.scatter(x[97, 165-128], y[97, 165-128], color='fuchsia', s=25)
    m.scatter(x[98, 185-128], y[98, 185-128], color='green', s=25)
    m.scatter(x[99, 205-128], y[99, 205-128], color='fuchsia', s=25)
    m.scatter(x[102, 232-128], y[102, 232-128], color='green', s=25)

    ax3 = fig.add_subplot(5, 1, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon_cam, meshlat_cam)
    m.contourf(x, y, c_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contourf(x, y, c_datas, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.85)
    m.scatter(x[73, 143-72], y[73, 143-72], color='fuchsia', s=25)
    m.scatter(x[72, 17+72], y[72, 17+72], color='green', s=25)
    m.scatter(x[72, 32+72], y[72, 32+72], color='fuchsia', s=25)
    m.scatter(x[71, 92-72], y[71, 92-72], color='fuchsia', s=25)
    m.scatter(x[71, 105-72], y[71, 105-72], color='green', s=25)
    m.scatter(x[71, 116-72], y[71, 116-72], color='fuchsia', s=25)
    m.scatter(x[71, 127-72], y[71, 127-72], color='green', s=25)

    ax4 = fig.add_subplot(5, 1, 4)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, g_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contourf(x, y, g_datas, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax4.set_yticks(np.arange(-30., 120., 30.))
    ax4.set_xticks(np.arange(-180., 270., 90.))
    ax4.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax4.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax4.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax4.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.85)
    m.scatter(x[30, 0+96], y[30, 0+96], color='fuchsia', s=25)
    m.scatter(x[32, 21+96], y[32, 21+96], color='green', s=25)
    m.scatter(x[31, 42+96], y[31, 42+96], color='fuchsia', s=25)
    m.scatter(x[37, 94-96], y[37, 94-96], color='green', s=25)
    m.scatter(x[39, 121-96], y[39, 121-96], color='fuchsia', s=25)
    m.scatter(x[34, 138-96], y[34, 138-96], color='green', s=25)
    m.scatter(x[35, 154-96], y[35, 154-96], color='fuchsia', s=25)
    m.scatter(x[30, 172-96], y[30, 172-96], color='green', s=25)

    ax5 = fig.add_subplot(5, 1, 5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, s_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                      vmax=np.max(ctrs), extend='both')
    m.contourf(x, y, s_datas, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax5.set_yticks(np.arange(-30., 120., 30.))
    ax5.set_xticks(np.arange(-180., 270., 90.))
    ax5.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax5.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax5.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax5.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax5.set_title('e', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.85)
    m.scatter(x[30, 0+96], y[30, 0+96], color='green', s=25)
    m.scatter(x[32, 21+96], y[32, 21+96], color='fuchsia', s=25)
    m.scatter(x[31, 42+96], y[31, 42+96], color='green', s=25)
    m.scatter(x[37, 94-96], y[37, 94-96], color='fuchsia', s=25)
    m.scatter(x[39, 121-96], y[39, 121-96], color='green', s=25)
    m.scatter(x[34, 138-96], y[34, 138-96], color='fuchsia', s=25)
    m.scatter(x[35, 154-96], y[35, 154-96], color='green', s=25)
    m.scatter(x[30, 172-96], y[30, 172-96], color='fuchsia', s=25)

    cbar_ax = fig.add_axes([0.4, 0.04, 0.2, 0.005])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Temperature ($^{\circ}$C)', size=12,
                fontsize=12, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    '''
    cbar_ax = fig.add_axes([0.552, 0.04, 0.42, 0.005])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=12,
                fontsize=12, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    '''
    pad = 30  # in points
    '''
    ax1.annotate('V200', xy=(0.5, 1), xytext=(0, 6), fontsize=16,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('Precipitation', xy=(0.5, 1), xytext=(0, 6), fontsize=16,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    '''
    ax1.annotate('HadAM3P', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax2.annotate('MIROC5', xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax2.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('CAM4', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax4.annotate('GHG', xy=(0, 0.5), xytext=(-ax4.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax4.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax5.annotate('SST', xy=(0, 0.5), xytext=(-ax5.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax5.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=0, top=.99, bottom=0.075,
                        left=0, right=1)


def fig4(had, mir, cam, ncep):
    from scipy import interpolate
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lat_cam = np.linspace(-90, 90, 96)
    lon_cam = np.linspace(0, 357.5, 144)
    lat_mir = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon_mir = np.arange(0, 360, 360/256)
    lon_ncep = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd.mon.mean.nc', 'lon')
    lat_ncep = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/NCEP-DOE/uwnd.mon.mean.nc', 'lat')
    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')

    h_data1 = had[5]#np.mean(had[5:8], axis=0)
    m_data1 = mir[5]#np.mean(mir[5:8], axis=0)
    c_data1 = cam[5]#np.mean(cam[5:8], axis=0)
    ncep1 = ncep[5]#np.mean(ncep[5:8], axis=0)
    
    h = interpolate.interp2d(lon, lat[::-1], h_data1)
    h_data = h(lon42, lat42[::-1])
    h = interpolate.interp2d(lon_mir, lat_mir[::-1], m_data1)
    m_data = h(lon42, lat42[::-1])
    h = interpolate.interp2d(lon_cam, lat_cam[::-1], c_data1)
    c_data = h(lon42, lat42[::-1])
    h = interpolate.interp2d(lon_ncep, lat_ncep[::-1], ncep1)
    n_data = h(lon42, lat42[::-1])

    h_data, lon1 = shiftgrid(180., h_data, lon42, start=False)
    h_data, lon1 = addcyclic(h_data, lon1)
    m_data, lon1 = shiftgrid(180., m_data, lon42, start=False)
    m_data, lon1 = addcyclic(m_data, lon1)
    c_data, lon1 = shiftgrid(180., c_data, lon42, start=False)
    c_data, lon1 = addcyclic(c_data, lon1)
    n_data, lon1 = shiftgrid(180., n_data, lon42, start=False)
    n_data, lon1 = addcyclic(n_data, lon1)

    meshlon, meshlat = np.meshgrid(lon1, lat42)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-16, 16, 17)
    ctrs2 = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1]) * 40
    ctrs2 = np.concatenate((np.linspace(-40, -5, 8), np.linspace(5, 40, 8)))
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = plt.subplot2grid((3, 6), (0, 0), colspan=5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, h_data-n_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contour(x, y, n_data, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax2 = plt.subplot2grid((3, 6), (1, 0), colspan=5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, m_data-n_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.contour(x, y, n_data, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax3 = plt.subplot2grid((3, 6), (2, 0), colspan=5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, c_data-n_data, ctrs, cmap=newcmap,
                      vmin=np.min(ctrs), vmax=np.max(ctrs), extend='both')
    m.contour(x, y, n_data, ctrs2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('e', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax4 = plt.subplot2grid((3, 6), (0, 5))
    for axis in ['top', 'bottom', 'left', 'right']:
        ax4.spines[axis].set_linewidth(2)
    ax4.plot(np.mean(h_data, axis=1), lat42, color='red', label='Model')
    ax4.plot(np.mean(n_data, axis=1), lat42, color='k', label='NCEP')
    ax4.plot(np.mean(h_data-n_data, axis=1), lat42, color='blue', label='Bias')
    ax4.set_ylim([-30, 90])
    ax4.set_xlim([-7.5, 7.5])
    ax4.set_yticks(np.arange(0., 90., 30.))
    ax4.set_xticks(np.arange(-10., 30., 5.))
    ax4.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax4.axvline(color='k', ls='--', lw=2)
    ax4.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.9)
    ax4.legend(loc=1)

    ax5 = plt.subplot2grid((3, 6), (1, 5))
    for axis in ['top', 'bottom', 'left', 'right']:
        ax5.spines[axis].set_linewidth(2)
    ax5.plot(np.mean(m_data, axis=1), lat42, color='red')
    ax5.plot(np.mean(n_data, axis=1), lat42, color='k')
    ax5.plot(np.mean(m_data-n_data, axis=1), lat42, color='blue')
    ax5.set_ylim([-30, 90])
    ax5.set_xlim([-7.5, 7.5])
    ax5.set_xticks(np.arange(-10., 30., 5.))
    ax5.set_yticks(np.arange(0., 90., 30.))
    ax5.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax5.axvline(color='k', ls='--', lw=2)
    ax5.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.9)

    ax6 = plt.subplot2grid((3, 6), (2, 5))
    for axis in ['top', 'bottom', 'left', 'right']:
        ax6.spines[axis].set_linewidth(2)
    ax6.plot(np.mean(c_data, axis=1), lat42, color='red')
    ax6.plot(np.mean(n_data, axis=1), lat42, color='k')
    ax6.plot(np.mean(c_data-n_data, axis=1), lat42, color='blue')
    ax6.set_ylim([-30, 90])
    ax6.set_xlim([-7.5, 7.5])
    ax6.set_xticks(np.arange(-5, 7.5, 2.5))
    ax6.set_xticks(np.arange(-10., 30., 5.))
    ax6.set_yticks(np.arange(0., 90., 30.))
    ax6.tick_params(labelleft='off', labelbottom='on', which='major',
                    direction='out', length=5, width=2)
    ax6.set_xlabel('Zonal mean u200 (ms$^{-1}$)', size=12)
    ax6.axvline(color='k', ls='--', lw=2)
    ax6.set_title('f', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.9)

    cbar_ax = fig.add_axes([0.233, 0.04, 0.4, 0.005])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Difference (ms$^{-1}$)', size=12,
                fontsize=12, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 37.5  # in points

    ax1.annotate('HadAM3P', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax2.annotate('MIROC5', xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax2.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('CAM4', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)

    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=-.2, top=.99, bottom=0.085,
                        left=0.06, right=0.95)


def fig5(f_ghg, v_ghg, f_sst, v_sst):
    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')

    g_data = (np.mean(v_ghg[16:21], axis=0)-v_ghg[0]) * 1e6
    s_data = (np.mean(v_sst[16:21], axis=0)-v_sst[0]) * 1e6
    g_dataf = f_ghg[5] * 1e11
    s_dataf = f_sst[5] * 1e11

    g_data, lon1 = shiftgrid(180., g_data, lon42, start=False)
    g_data, lon1 = addcyclic(g_data, lon1)
    s_data, lon1 = shiftgrid(180., s_data, lon42, start=False)
    s_data, lon1 = addcyclic(s_data, lon1)
    g_dataf, lon1 = shiftgrid(180., g_dataf, lon42, start=False)
    g_dataf, lon1 = addcyclic(g_dataf, lon1)
    s_dataf, lon1 = shiftgrid(180., s_dataf, lon42, start=False)
    s_dataf, lon2 = addcyclic(s_dataf, lon1)

    meshlon42, meshlat42 = np.meshgrid(lon2, lat42)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-1, 1, 17)
    ctrs2 = np.linspace(-2, 2, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(2, 2, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=30,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon42, meshlat42)
    m.contourf(x, y, g_dataf, ctrs2, cmap=newcmap, vmin=np.min(ctrs2),
               vmax=np.max(ctrs2), extend='both')
    ax1.axhline(15, color='Fuchsia', lw=2)
    ax1.axhline(0, color='Fuchsia', lw=2)
    m.drawparallels(np.arange(-30, 45, 15),
                    labels=[True, False, False, True], linewidth=0, xoffset=10,
                    fontsize=12)
    ax1.set_yticks(np.arange(-30., 45., 15.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 35., 5.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax2 = fig.add_subplot(2, 2, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=30,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon42, meshlat42)
    plot = m.contourf(x, y, s_dataf, ctrs2, cmap=newcmap, vmin=np.min(ctrs2),
                      vmax=np.max(ctrs2), extend='both')
    ax2.axhline(15, color='Fuchsia', lw=2)
    ax2.axhline(0, color='Fuchsia', lw=2)
    ax2.set_yticks(np.arange(-30., 45., 15.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 35., 5.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax3 = fig.add_subplot(2, 2, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon42, meshlat42)
    m.contourf(x, y, g_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10,
                    fontsize=12)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10,
                    fontsize=12)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    ax4 = fig.add_subplot(2, 2, 4)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon42, meshlat42)
    plot1 = m.contourf(x, y, s_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                       vmax=np.max(ctrs), extend='both')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10,
                    fontsize=12)
    ax4.set_yticks(np.arange(-30., 120., 30.))
    ax4.set_xticks(np.arange(-180., 270., 90.))
    ax4.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax4.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax4.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax4.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=24,
                  fontweight='bold', y=0.99)

    cbar_ax = fig.add_axes([0.3, 0.14, 0.4, 0.0075])
    b = fig.colorbar(plot1, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Forcing (2x10$^{-11}$ s$^{-2}$) & meridional wind response (ms$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16) 
    ax1.annotate('GHG', xy=(0.5, 1), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('SST', xy=(0.5, 1), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')

    plt.show()
    plt.subplots_adjust(hspace=-.5, wspace=0.05, top=.91, bottom=0.1,
                        left=0.05, right=0.95)


def fig6(w_ghg_j, w_ghg_j_ncep, w_sst_j, w_sst_j_ncep,
          v_greens, v_greens_ncep, f_ghg, f_sst):
    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')
    lons = np.arange(-180, 180+5.625, 360/64)

    #f_ghg_m = np.mean(f_ghg[4, 28:32, :], axis=0)
    #f_sst_m = np.mean(f_sst[4, 28:32, :], axis=0)
    f_ghg_m = np.mean(f_ghg[5, 26:29, :].transpose()*np.cos(lat42[26:29]*np.pi/180), axis=1)/np.mean(np.cos(lat42[26:29]*np.pi/180))
    f_sst_m = np.mean(f_sst[5, 26:29, :].transpose()*np.cos(lat42[26:29]*np.pi/180), axis=1)/np.mean(np.cos(lat42[26:29]*np.pi/180))
    f_ghg_m = np.concatenate((f_ghg_m[64:], f_ghg_m[:65]))
    f_sst_m = np.concatenate((f_sst_m[64:], f_sst_m[:65]))*-1

    v_datai = v_greens[10]
    v_datap = v_greens[61]*-1  # or2??
    v_dataa = v_greens[56]
    v_dataall = v_datai + v_dataa + v_datap

    v_datain = v_greens_ncep[10]
    v_datapn = v_greens_ncep[61]*-1  # or2??
    v_dataan = v_greens_ncep[56]
    v_dataalln = v_datain + v_dataan + v_datapn

    v_datai, lon1 = shiftgrid(180., v_datai, lon42, start=False)
    v_datai, lon1 = addcyclic(v_datai, lon1)
    v_datap, lon1 = shiftgrid(180., v_datap, lon42, start=False)
    v_datap, lon1 = addcyclic(v_datap, lon1)
    v_dataa, lon1 = shiftgrid(180., v_dataa, lon42, start=False)
    v_dataa, lon1 = addcyclic(v_dataa, lon1)
    v_dataall, lon1 = shiftgrid(180., v_dataall, lon42, start=False)
    v_dataall, lon1 = addcyclic(v_dataall, lon1)

    v_datain, lon1 = shiftgrid(180., v_datain, lon42, start=False)
    v_datain, lon1 = addcyclic(v_datain, lon1)
    v_datapn, lon1 = shiftgrid(180., v_datapn, lon42, start=False)
    v_datapn, lon1 = addcyclic(v_datapn, lon1)
    v_dataan, lon1 = shiftgrid(180., v_dataan, lon42, start=False)
    v_dataan, lon1 = addcyclic(v_dataan, lon1)
    v_dataalln, lon1 = shiftgrid(180., v_dataalln, lon42, start=False)
    v_dataalln, lon1 = addcyclic(v_dataalln, lon1)

    meshlon, meshlat = np.meshgrid(lon1, lat42)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    ax1 = plt.subplot2grid((6, 2), (0, 0))
    ax2 = plt.subplot2grid((6, 2), (1, 0))
    ax3 = plt.subplot2grid((6, 2), (0, 1), rowspan=2)
    ax4 = plt.subplot2grid((6, 2), (2, 0), rowspan=2)
    ax5 = plt.subplot2grid((6, 2), (2, 1), rowspan=2)
    ax6 = plt.subplot2grid((6, 2), (4, 0), rowspan=2)
    ax7 = plt.subplot2grid((6, 2), (4, 1), rowspan=2)

    ctrs = np.linspace(-1, 1, 17)*.4
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    #ax1 = fig.add_subplot(3, 2, 1)
    w_ghg_j = np.concatenate((w_ghg_j[32:], w_ghg_j[:33]))
    w_ghg_j_ncep = np.concatenate((w_ghg_j_ncep[32:], w_ghg_j_ncep[:33]))
    # w_ghg_jja = np.concatenate((w_ghg_jja[32:], w_ghg_jja[:33]))
    ax1.plot(lons, w_ghg_j, color='#1258DC', label='GHG')
    ax1.plot(lons, w_ghg_j_ncep, color='#7FBD32', label='NCEP')
    # ax1.plot(lons, w_ghg_a, color='#1258DC', label='Aug')
    # ax1.plot(lons, w_ghg_jja, color='#091834', label='JJA')
    ax1.plot(np.arange(-180, 180+2.8125, 2.8125),
             f_ghg_m*np.max(w_ghg_j)/np.max(f_ghg_m),
             color='#FD4D0C', label='Forcing')
    ax1.axhline(0, color='k')
    ax1.axhline(-.14, color='k', lw=3)
    ax1.axhline(.14, color='k', lw=4)
    ax1.axvline(-180, color='k', lw=4)
    ax1.axvline(180, color='k', lw=3)
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_yticks(np.arange(-0.1, .20, .05))
    ax1.tick_params(length=5, width=2)
    ax1.tick_params(labelsize=12) 
    ax1.set_xlim([-180, 180])
    ax1.set_ylim([-0.14, .14])
    ax1.legend(loc=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.98, x=0.005)

    #ax2 = fig.add_subplot(3, 2, 2)
    w_sst_j = np.concatenate((w_sst_j[32:], w_sst_j[:33]))*-1
    w_sst_j_ncep = np.concatenate((w_sst_j_ncep[32:], w_sst_j_ncep[:33]))*-1
    # w_sst_jja = np.concatenate((w_sst_jja[32:], w_sst_jja[:33]))*-1
    ax2.plot(lons, w_sst_j, color='#1258DC', label='SST')
    ax2.plot(lons, w_sst_j_ncep, color='#7FBD32', label='NCEP')
    # ax2.plot(lons, w_sst_a, color='#1258DC', label='Aug')
    # ax2.plot(lons, w_sst_jja, color='#091834', label='JJA')
    ax2.plot(np.arange(-180, 180+2.8125, 2.8125),
             f_sst_m*np.max(w_sst_j)/np.max(f_sst_m),
             color='#FD4D0C', label='Forcing')
    ax2.axhline(0, color='k')
    ax2.axhline(-.05, color='k', lw=3)
    ax2.axhline(.05, color='k', lw=4)
    ax2.axvline(-180, color='k', lw=4)
    ax2.axvline(180, color='k', lw=3)
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_yticks(np.arange(-0.04, .08, .02))
    ax2.tick_params(labelsize=12) 
    ax2.set_xlim([-180, 180])
    ax2.set_ylim([-0.05, .05])
    ax2.legend(loc=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.97, x=0.005)

    #ax3 = fig.add_subplot(3, 2, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax3)
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, v_dataa, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    ax3.scatter(lons[56]-180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.98, x=0.005)

    #ax4 = fig.add_subplot(3, 2, 4)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax4)
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, v_datap, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    ax4.scatter(lons[61]-180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10,
                    fontsize=12)
    ax4.set_yticks(np.arange(-30., 120., 30.))
    ax4.set_xticks(np.arange(-180., 270., 90.))
    ax4.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax4.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax4.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax4.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.98, x=0.005)

    #ax5 = fig.add_subplot(3, 2, 5)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax5)
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, v_datai, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                      vmax=np.max(ctrs), extend='both')
    ax5.scatter(lons[10]+180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    ax5.set_yticks(np.arange(-30., 120., 30.))
    ax5.set_xticks(np.arange(-180., 270., 90.))
    ax5.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax5.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax5.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax5.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax5.set_title('e', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.98, x=0.005)

    #ax6 = fig.add_subplot(3, 2, 6)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax6)
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, v_dataall, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    ax6.scatter(lons[10]+180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    ax6.scatter(lons[56]-180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    ax6.scatter(lons[61]-180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10,
                    fontsize=12)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10,
                    fontsize=12)
    ax6.set_yticks(np.arange(-30., 120., 30.))
    ax6.set_xticks(np.arange(-180., 270., 90.))
    ax6.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax6.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax6.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax6.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax6.set_title('f', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=0.98, x=0.005)

    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax7)
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    m.contourf(x, y, v_dataalln, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    ax7.scatter(lons[10]+180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    ax7.scatter(lons[56]-180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    ax7.scatter(lons[61]-180, lat42[29], marker='x',
                color='Fuchsia', linewidth=2, s=50)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10,
                    fontsize=12)
    ax7.set_yticks(np.arange(-30., 120., 30.))
    ax7.set_xticks(np.arange(-180., 270., 90.))
    ax7.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax7.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax7.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax7.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax7.set_title('g', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', y=1, x=0.005)

    cbar_ax = fig.add_axes([0.94, 0.3, 0.0075, 0.4])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Meridional wind response (ms$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)

    ax3.annotate('Western Atlantic', xy=(0.5, .98), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax4.annotate('Eastern Atlantic', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax5.annotate('Indian Ocean', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax6.annotate('Combined response using NAT background', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax7.annotate('Combined response using NCEP background', xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')

    plt.show()
    plt.subplots_adjust(hspace=0.4, wspace=0.05, top=.97, bottom=0.04,
                        left=0.05, right=0.925)


def postage(va, pr):
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    ctrs = np.linspace(-1, 1, 17)
    ctrs2 = np.array([-2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6])
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap1 = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap1)

    fig, axs = plt.subplots(3, 3, facecolor='w', edgecolor='k', linewidth=2)

    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    names = ['CNRM-CM5', 'CanAM4', 'FGOALS-g2', 'HadGEM2-A', 'IPSL-CM5A-LR', 'IPSL-CM5B-LR', 'MPI-ESM-LR', 'MPI-ESM-MR', 'bcc-csm1-1']
    for i in range(9):
        v_data, lon1 = shiftgrid(180., va[i], lon, start=False)
        v_data, lon1 = addcyclic(v_data, lon1)
        p_data, lon1 = shiftgrid(180., pr[i], lon, start=False)
        p_data, lon1 = addcyclic(p_data, lon1)
        meshlon, meshlat = np.meshgrid(lon1, lat)
        ax1 = axs[int(np.floor(i/3)), np.remainder(i, 3)]
        #for axis in ['top', 'bottom', 'left', 'right']:
         #   ax1.spines[axis].set_linewidth(2)
        m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                    llcrnrlon=-180, urcrnrlon=180, resolution='c', ax=ax1)
        m.drawcoastlines(color='gray')
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        plot = m.contourf(x, y, p_data, ctrs*4, cmap=newcmap, vmin=np.min(ctrs*4),
                          vmax=np.max(ctrs*4), extend='both')
        m.contour(x, y, v_data, ctrs2*4, colors='k')
        if np.remainder(i, 3) == 0:
            m.drawparallels(np.arange(0, 90, 30),
                            labels=[True, False, False, True], linewidth=0, xoffset=10)
        if i > 5:
            m.drawmeridians(np.arange(-180, 270, 90),
                            labels=[True, False, False, True], linewidth=0, yoffset=10)
        ax1.set_yticks(np.arange(-30., 120., 30.))
        ax1.set_xticks(np.arange(-180., 270., 90.))
        ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                        direction='out', length=5, width=2)
        ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
        ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
        ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                        direction='out', length=4, width=1)
        ax1.annotate(names[i], xy=(0.5, .995), xytext=(0, 6), fontsize=24,
                     fontname='Arial',
                     xycoords='axes fraction', textcoords='offset points',
                     ha='center', va='baseline')
    cbar_ax = fig.add_axes([0.94, 0.3, 0.0075, 0.4])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=0)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    plt.subplots_adjust(hspace=0, wspace=0.08, top=.97, bottom=0.03,
                        left=0.03, right=.925)


def fig8(sst, ghg, amipco2, amip4K, v_sst, v_ghg,
         sstp, ghgp, amipco2p, amip4Kp, f_sst, f_ghg):
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')

    b_data = ((np.mean(v_ghg[16:21], axis=0)-v_ghg[0]) -
              (np.mean(v_sst[16:21], axis=0)-v_sst[0])) * 1e6
    b_dataf = (f_ghg[5] - f_sst[5]) * 1e11

    b_data, lon2 = shiftgrid(180., b_data, lon42, start=False)
    b_data, lon2 = addcyclic(b_data, lon2)
    b_dataf, lon2 = shiftgrid(180., b_dataf, lon42, start=False)
    b_dataf, lon2 = addcyclic(b_dataf, lon2)

    g_data = np.mean(ghg[5:8], axis=0) - np.mean(sst[5:8], axis=0)
    g_datap = np.mean(ghgp[5:8], axis=0) - np.mean(sstp[5:8], axis=0)

    g_data, lon1 = shiftgrid(180., g_data, lon, start=False)
    g_data, lon1 = addcyclic(g_data, lon1)
    g_datap, lon1 = shiftgrid(180., g_datap, lon, start=False)
    g_datap, lon1 = addcyclic(g_datap, lon1)

    a_data, lon1 = shiftgrid(180., amipco2.mean(axis=0)-amip4K.mean(axis=0),
                             lon, start=False)
    a_data, lon1 = addcyclic(a_data, lon1)

    ap_data, lon1 = shiftgrid(180., amipco2p.mean(axis=0)-amip4Kp.mean(axis=0),
                              lon, start=False)
    ap_data, lon1 = addcyclic(ap_data, lon1)

    meshlon, meshlat = np.meshgrid(lon1, lat)
    meshlon42, meshlat42 = np.meshgrid(lon2, lat42)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-1, 1, 17)
    ctrs2 = np.array([-2.6, -2.4, -2.2, -2.0, -1.8, -1.6, -1.4, -1.2, -1, -.8,
                      -.6, -.4, -.2, .2, .4, .6, .8, 1, 1.2, 1.4, 1.6, 1.8, 2,
                      2.2, 2.4, 2.6])

    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap1 = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    my_cmap = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(3, 1, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot1 = m.contourf(x, y, ap_data, ctrs*6, cmap=newcmap, vmin=np.min(ctrs*6),
                       vmax=np.max(ctrs*6), extend='both')
    m.contour(x, y, a_data, ctrs2*4, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax2 = fig.add_subplot(3, 1, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot2 = m.contourf(x, y, g_datap, ctrs*2, cmap=newcmap, vmin=np.min(ctrs*2),
                       vmax=np.max(ctrs*2), extend='both')
    m.contour(x, y, g_data, ctrs2*2, colors='k')
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax3 = fig.add_subplot(3, 1, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon42, meshlat42)
    plot3 = m.contourf(x[26:32], y[26:32], b_dataf[26:32], ctrs*3, cmap=newcmap1, vmin=np.min(ctrs*3),
                       vmax=np.max(ctrs*3), extend='both')
    m.contour(x, y, b_data, ctrs2*1.5, colors='k')
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    cbar_ax = fig.add_axes([0.8, 0.7125, 0.0075, 0.25])
    b = fig.colorbar(plot1, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    cbar_ax = fig.add_axes([0.8, 0.385, 0.0075, 0.25])
    b = fig.colorbar(plot2, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    cbar_ax = fig.add_axes([0.8, 0.0575, 0.0075, 0.25])
    b = fig.colorbar(plot3, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label='Difference (1x10$^{-11}$ s$^{-2}$)', size=16,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'ymajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)
    cbar_ax.tick_params(labelsize=16)
    pad = 40  # in points

    ax1.annotate('CMIP5 AMIP', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax2.annotate('HadAM3P', xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax2.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('Barotropic', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)

    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=0, top=.99, bottom=0.04,
                        left=0, right=1)


def fig_rean(data, sig, gis, sigt, lat, lon):

    latg = np.linspace(-89, 89, 90)
    long = np.linspace(-179, 179, 180)
    gis, long1 = addcyclic(gis, long)
    s_data, long = addcyclic(sigt, long)
    meshlon1, meshlat1 = np.meshgrid(long, latg)

    h_data, lon1 = shiftgrid(180., data, lon, start=False)
    h_data, lon1 = addcyclic(h_data, lon1)
    m_data, lon1 = shiftgrid(180., sig, lon, start=False)
    m_data, lon1 = addcyclic(m_data, lon1)
    meshlon, meshlat = np.meshgrid(lon1, lat)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-.2, .2, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(2, 1, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, h_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
                      vmax=np.max(ctrs), extend='both')
    m.contourf(x, y, m_data, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    #m.drawmeridians(np.arange(-180, 270, 90),
     #               labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('NCEP2 1979-2015', fontsize=20, fontname='Arial')
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax2 = fig.add_subplot(2, 1, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon1, meshlat1)
    plot1 = m.contourf(x, y, gis, ctrs/2.5, cmap=newcmap, vmin=np.min(ctrs/2.5),
                       vmax=np.max(ctrs/2.5), extend='both')
    m.contourf(x, y, s_data, levels=[0, .5, 1.5], hatches=["", "///"],
               alpha=0)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('GISTEMP 1979-2015', fontsize=20, fontname='Arial')
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    cbar_ax = fig.add_axes([0.875, 0.561, 0.01, 0.3])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.1f')
    b.set_label(label=r'v200 trend (ms$^{-1}$ yr$^{-1}$)', size=12,
                fontsize=12, fontname='Arial', labelpad=4)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    cbar_ax = fig.add_axes([0.875, 0.142, 0.01, 0.3])
    b = fig.colorbar(plot1, cax=cbar_ax, spacing='proportional',
                     orientation='vertical', extend='max', format='%.2f')
    b.set_label(label=r'Surface temperature trend (K yr$^{-1}$)', size=12,
                fontsize=12, fontname='Arial', labelpad=-1)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.yaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    plt.show()
    plt.subplots_adjust(hspace=0.1, wspace=0, top=.9, bottom=0.1,
                        left=0.1, right=.9)


def fig1su(had, mir, cam):
    def ttest(series1, series2, siglevel=5, testtype='two'):
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
        sigarray = np.full(np.shape(series1[0]), siglevel)
        if testtype == 'two':
            a = 2
        elif testtype == 'one':
            a = 1
        else:
            print("Error, test type must be 'one' or 'two'")
        d = np.sqrt(np.var(series1, axis=0, ddof=1)/len(series1) + np.var(series2, axis=0,
                                                             ddof=1)/len(series2))
        z1 = (np.mean(series1, axis=0) - np.mean(series2, axis=0)) / d
        p1 = 1 - st.norm.cdf(np.abs(z1))
        sig = np.greater_equal(sigarray, p1*100*a).astype(int)
        return sig

    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    lat_cam = np.linspace(-90, 90, 96)
    lon_cam = np.linspace(0, 357.5, 144)
    lat_mir = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    lon_mir = np.arange(0, 360, 360/256)
    h_data = (had['Plus15-Future_HCO2'][1].mean(axis=0)-had['Plus15-Future_LCO2'][1].mean(axis=0))*86400
    m_data = (mir['Plus15-Future_HCO2'][1].mean(axis=0)-mir['Plus15-Future_LCO2'][1].mean(axis=0))*86400
    c_data = (cam['Plus15-Future'][1].mean(axis=0)-cam['Plus15-Future_LCO2'][1].mean(axis=0))*2.5*86400

    h_data, lon1 = shiftgrid(180., h_data, lon, start=False)
    h_data, lon1 = addcyclic(h_data, lon1)
    m_data, lon1_mir = shiftgrid(180., m_data, lon_mir, start=False)
    m_data, lon1_mir = addcyclic(m_data, lon1_mir)
    c_data, lon1_cam = shiftgrid(180., c_data, lon_cam, start=False)
    c_data, lon1_cam = addcyclic(c_data, lon1_cam)

    h_sig = ttest(had['Plus15-Future_HCO2'][1], had['Plus15-Future_LCO2'][1])
    m_sig = ttest(mir['Plus15-Future_HCO2'][1], mir['Plus15-Future_LCO2'][1])
    c_sig = ttest(cam['Plus15-Future'][1], cam['Plus15-Future_LCO2'][1])

    h_sig, lon1 = shiftgrid(180., h_sig, lon, start=False)
    h_sig, lon1 = addcyclic(h_sig, lon1)
    m_sig, lon1_mir = shiftgrid(180., m_sig, lon_mir, start=False)
    m_sig, lon1_mir = addcyclic(m_sig, lon1_mir)
    c_sig, lon1_cam = shiftgrid(180., c_sig, lon_cam, start=False)
    c_sig, lon1_cam = addcyclic(c_sig, lon1_cam)

    meshlon, meshlat = np.meshgrid(lon1, lat)
    meshlon_mir, meshlat_mir = np.meshgrid(lon1_mir, lat_mir)
    meshlon_cam, meshlat_cam = np.meshgrid(lon1_cam, lat_cam)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ctrs = np.linspace(-1, 1, 17)
    #ctrs2 = np.array([-1, -.875, -.75, -.625, -.5, -.375, -.25, -.125, .125, .25, .375, .5, .625, .75, .875, 1])
    ctrs2 = np.array([-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1])
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    my_cmap = my_cmap[::-1]
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    ax1 = fig.add_subplot(3, 1, 1)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, h_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    plt.contourf(x, y, h_sig, levels=[0, .5, 1.5],
                 hatches=["", ".."], alpha=0)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax1.set_yticks(np.arange(-30., 120., 30.))
    ax1.set_xticks(np.arange(-180., 270., 90.))
    ax1.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax1.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax1.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax1.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax2 = fig.add_subplot(3, 1, 2)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon_mir, meshlat_mir)
    m.contourf(x, y, m_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    plt.contourf(x, y, m_sig, levels=[0, .5, 1.5],
                 hatches=["", ".."], alpha=0)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax2.set_yticks(np.arange(-30., 120., 30.))
    ax2.set_xticks(np.arange(-180., 270., 90.))
    ax2.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax2.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax2.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax2.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    ax3 = fig.add_subplot(3, 1, 3)
    m = Basemap(projection='cyl', llcrnrlat=-30, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines(color='gray')
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon_cam, meshlat_cam)
    m.contourf(x, y, c_data, ctrs, cmap=newcmap, vmin=np.min(ctrs),
               vmax=np.max(ctrs), extend='both')
    plt.contourf(x, y, c_sig, levels=[0, .5, 1.5],
                 hatches=["", ".."], alpha=0)
    m.drawmeridians(np.arange(-180, 270, 90),
                    labels=[True, False, False, True], linewidth=0, yoffset=10)
    m.drawparallels(np.arange(0, 90, 30),
                    labels=[True, False, False, True], linewidth=0, xoffset=10)
    ax3.set_yticks(np.arange(-30., 120., 30.))
    ax3.set_xticks(np.arange(-180., 270., 90.))
    ax3.tick_params(labelleft='off', labelbottom='off', which='major',
                    direction='out', length=5, width=2)
    ax3.set_yticks(np.arange(-30., 105., 15.), minor=True)
    ax3.set_xticks(np.arange(-180., 210., 30.), minor=True)
    ax3.tick_params(labelleft='off', labelbottom='off', which='minor',
                    direction='out', length=4, width=1)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.04, y=0.9)

    cbar_ax = fig.add_axes([0.4, 0.04, 0.2, 0.005])
    b = fig.colorbar(plot, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='Difference (mm day$^{-1}$)', size=12,
    #b.set_label(label='Difference (ms$^{-1}$)', size=12,
                fontsize=12, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=12)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    pad = 40  # in points

    ax1.annotate('HadAM3P', xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax2.annotate('MIROC5', xy=(0, 0.5), xytext=(-ax2.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax2.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)
    ax3.annotate('CAM4', xy=(0, 0.5), xytext=(-ax3.yaxis.labelpad - pad, 0),
                 fontsize=16, fontname='Arial',
                 xycoords=ax3.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)

    plt.show()
    plt.subplots_adjust(hspace=0.075, wspace=0, top=.99, bottom=0.085,
                        left=0, right=1)



















