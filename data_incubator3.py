#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 16:54:14 2018

@author: bakerh
"""

def ncread(filelocation, invariable):
    '''
    ncread outputs numpy arrays of called variables in netCDF4 file

    Parameters
    ----------
    fileloaction : str
        NetCDFfile
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variable : array
        Array of variable specified
    '''
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    nc_f = filelocation  # filename
    nc_fid = Dataset(nc_f, 'r')
    variable = nc_fid.variables[invariable][:]
    return variable


def main1():
    from scipy import interpolate
    import glob
    import numpy as np



    def powerfunction(w100):
        power = np.zeros((np.shape(w100)[0], np.shape(w100)[1]))
        # loop over grid
        for i in range(np.shape(w100)[0]):
            for j in range(np.shape(w100)[1]):
                if w100[i, j] < 2.5:
                    power[i, j] = 0
                elif w100[i, j] <= 12 and w100[i, j] >= 2.5:
                    # power curve from IEC 61400-12-1. Wind Turbines - Part 12-1
                    power[i, j] = .75*(-.05*w100[i, j]**5+1.24*w100[i, j]**4-9.74*w100[i, j]**3+45.32*w100[i, j]**2-78.08*w100[i, j]+35.62)
                elif w100[i, j] >=12 and w100[i, j] <= 25.5:
                    power[i, j] = 1500
                else:
                    power[i, j] = 0
        return power

    def wind_100(w, t, p):
        # compute 100m wind from 10m wind
        R = 287.058  # specific gas constant of dry air
        w100 = w*10**(1/7)  # wind scaling power law with altitude doi:10.1115/1.3035818
        rho = p/(R*t)  # compute density
        w100_s = w100*(rho/1.225)**(1/3)  # scale wind for density IEC 61400-12-1
        return w100_s

    def shiftlat(grid, lat, lon, latw):
        # regrid wind to same lat lon grid as temperature
        regrid = np.zeros((len(lat), len(lon)))
        for i in range(np.ma.size(latw)-1):
            regrid[i+1, :] = ((grid[i, :]*np.cos(latw[i]*np.pi/180) +
                                  grid[i+1, :]*np.cos(latw[i+1]*np.pi/180)) /
                                 (2*np.cos(lat[i+1]*np.pi/180)))
        return regrid

    # create lat lon coords to interpolate individual model data on to
    lat_n96 = np.arange(90, -91.25, -1.25)
    lon_n96 = np.arange(0, 360, 1.875)
    # model list
    models = ['ACCESS1-0_', 'ACCESS1-3_', 'BNU-ESM_', 'CanESM2_',
              'CESM1-CAM5_', 'CMCC-CESM_', 'CMCC-CMS_', 'CNRM-CM5_',
              'CSIRO-Mk3-6-0_', 'GFDL-CM3_', 'GFDL-ESM2G_',
              'GFDL-ESM2M_', 'GISS-E2-H-CC_', 'GISS-E2-R-CC_',
              'GISS-E2-R_', 'HadGEM2-AO', 'HadGEM2-CC', 'HadGEM2-ES',
              'IPSL-CM5A-LR_', 'IPSL-CM5A-MR_',
              'IPSL-CM5B-LR_', 'MPI-ESM-LR_', 'MPI-ESM-MR_']

    p_mmm_hist = np.zeros((145, 192))
    p_mmm_rcp85 = np.zeros((145, 192))
    # loop over models calculating wind power in each model for historical and future scenarios
    for j, model in enumerate(models):
        # load variables for particular model
        a1 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/historical/tas/tas_Amon_' + model + '*')
        a2 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/historical/psl/psl_Amon_' + model + '*')
        a3 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/historical/sfcWind/sfcWind_Amon_' + model + '*')
        b1 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/rcp85/tas/tas_Amon_' + model + '*')
        b2 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/rcp85/psl/psl_Amon_' + model + '*')
        b3 = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/CMIP5/rcp85/sfcWind/sfcWind_Amon_' + model + '*')
        tas_hist = ncread(a1[0], 'tas')
        psl_hist = ncread(a2[0], 'psl')
        sfcWind_hist = ncread(a3[0], 'sfcWind')
        tas_rcp85 = ncread(b1[0], 'tas')
        psl_rcp85 = ncread(b2[0], 'psl')
        sfcWind_rcp85 = ncread(b3[0], 'sfcWind')
        # some models are split across multiple files, combine here
        if len(a1) > 1:
            for i in range(len(a1)-1):
                tas_hist = np.concatenate((tas_hist, ncread(a1[i+1], 'tas')), axis=0)
        if len(a2) > 1:
            for i in range(len(a2)-1):
                psl_hist = np.concatenate((psl_hist, ncread(a2[i+1], 'psl')), axis=0)
        if len(a3) > 1:
            for i in range(len(a3)-1):
                sfcWind_hist = np.concatenate((sfcWind_hist, ncread(a3[i+1], 'sfcWind')), axis=0)
        if len(b1) > 1:
            for i in range(len(b1)-1):
                tas_rcp85 = np.concatenate((tas_rcp85, ncread(b1[i+1], 'tas')), axis=0)
        if len(a2) > 1:
            for i in range(len(b2)-1):
                psl_rcp85 = np.concatenate((psl_rcp85, ncread(b2[i+1], 'psl')), axis=0)
        if len(b3) > 1:
            for i in range(len(b3)-1):
                sfcWind_rcp85 = np.concatenate((sfcWind_rcp85, ncread(b3[i+1], 'sfcWind')), axis=0)
        # pick 1980-2000 for hist runs and 2080-2100 for future scenario
        tas_hist = np.squeeze(np.mean(tas_hist[1560:1800], axis=0))
        psl_hist = np.squeeze(np.mean(psl_hist[1560:1800], axis=0))
        sfcWind_hist = np.squeeze(np.mean(sfcWind_hist[1560:1800], axis=0))
        tas_rcp85 = np.squeeze(np.mean(tas_rcp85[888:1128], axis=0))
        psl_rcp85 = np.squeeze(np.mean(psl_rcp85[888:1128], axis=0))
        sfcWind_rcp85 = np.squeeze(np.mean(sfcWind_rcp85[888:1128], axis=0))
        lat = ncread(a1[0], 'lat')
        lon = ncread(a1[0], 'lon')
        latw = ncread(a3[0], 'lat')
        # regrid if wind grid is different
        if len(lat) != len(latw):
            sfcWind_hist = shiftlat(sfcWind_hist, lat, lon, latw)
            sfcWind_rcp85 = shiftlat(sfcWind_rcp85, lat, lon, latw)

        # compute wind power in hist and future
        w100_hist = wind_100(sfcWind_hist, tas_hist, psl_hist)
        p_hist = powerfunction(w100_hist)
        w100_rcp85 = wind_100(sfcWind_rcp85, tas_rcp85, psl_rcp85)
        p_rcp85 = powerfunction(w100_rcp85)
        # interpolate to common grid
        f = interpolate.interp2d(lon, lat[::-1], p_hist)
        g = interpolate.interp2d(lon, lat[::-1], p_rcp85)
        p_hist_n96 = f(lon_n96, lat_n96)
        p_rcp85_n96 = g(lon_n96, lat_n96)

        p_mmm_hist += p_hist_n96
        p_mmm_rcp85 += p_rcp85_n96
        print(model)
    p_mmm_hist /= len(models)
    p_mmm_rcp85 /= len(models)
    return p_mmm_hist, p_mmm_rcp85


def maplot(pdata, colormax=1, colormin=-999, title='', anom='no'):
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
    # land-sea mask
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    pdata = np.ma.masked_array(pdata, np.logical_not(lsm))
    lat = np.arange(90, -91.25, -1.25)
    lon = np.arange(0, 360, 1.875)
    if colormin == -999:
        colormin = -colormax
    plt.figure()
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
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    cmap1 = newcmap
    if anom == 'no':
        cmap1 = 'plasma_r'
    ctrs = np.linspace(colormin, colormax, 17)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=cmap1, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    b = plt.colorbar(plot, orientation='horizontal',
                     aspect=50, shrink=0.75, spacing='proportional')
    b.set_label(label='Power (kW)')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()

p_mmm_hist, p_mmm_rcp85 = main1()
maplot(p_mmm_hist, 400, 0, title='Simulated wind power climatology')
maplot(p_mmm_rcp85-p_mmm_hist, 50, anom='yes', title='Predicted changes in wind power')
