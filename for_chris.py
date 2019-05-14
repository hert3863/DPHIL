#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 13:59:09 2019

@author: bakerh
"""

'''
README
All the functions require the inputs to have longitude -180 to 180 and latitude
90 to -90. They all output winter and summer in the same array
(first dimension is the season with winter first).

The first function computes the first EOF of MSLP, and you'll need to change
line 41 in the code to point to your control run directory. You'll also need
to change line 47 to have the correct variable name for MSLP in
your NetCDF files.

The second function computes the GTO. You'll need to change lines
89 and 90 to point at your control and forced run directories respectively,
and lines 98 and 100 to have the correct variable names again.

The final function will allow you to reconstruct the NAO using different SST
basins. It requires the observed NAO computed from the NOAA20thC reanalysis
as an input (attached as a csv). For this function you'll need the HadISST
dataset and you'll need to change line 150 to point to this.
'''


import numpy as np


def compute_eof(lat, lon, n_eof=1):
    # this function computes the first EOF of the model control runs
    # over region 20N-80N, 70W-40E
    # need data to be lon -180 to 180, lat 90 to -90
    # pass in n_eof=1 for first eof, n_eof=2 for second etc.
    import glob
    # point this towards control run directory
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/*')
    mslp_w = np.zeros((len(a), len(lat), len(lon)))
    mslp_s = np.zeros((len(a), len(lat), len(lon)))
    for i, item in enumerate(a):
        # change 'item16222_monthly_mean' to MSLP pressure
        # variable name in netcdf files
        mslp_96 = ncread(item, 'item16222_monthly_mean')[:, 0, :] / 100
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        # pulls out seasonal mean fields
        mslp_w[i] = np.mean(mslp_96[windices], axis=0)
        mslp_s[i] = np.mean(mslp_96[sindices], axis=0)
        print(str(i))
    eofs = np.zeros((2, len(lat), len(lon)))
    var = np.zeros((2, 3))
    # get EOF bounding indices
    r1 = np.abs(lat-80).argmin()
    r2 = np.abs(lat-20).argmin() + 1
    r3 = np.abs(lon+70).argmin()
    r4 = np.abs(lon-40).argmin() + 1
    # compute eofs
    eofs[0], var[0] = eof_clim(mslp_w, lat, lon, n_eof,
                               region=[r1, r2, r3, r4])
    eofs[1], var[1] = eof_clim(mslp_s, lat, lon, n_eof,
                               region=[r1, r2, r3, r4])
    return eofs


def compute_gto(eof, lat, lon):
    # function computes NAO GTO and outputs two fields,
    # the GTO and the GTO with significant values only
    # data must be organised with lon coords from -180 to 180
    # and lat coords from 90 to -90
    def win_sum(field):
        # this function computes the seasonl MSLP fields
        sindices = [17, 18, 19]
        windices = [11, 12, 13]
        field_w = np.mean(field[windices], axis=0)
        field_s = np.mean(field[sindices], axis=0)
        return field_w, field_s
    # get EOF bounding indices
    r1 = np.abs(lat-80).argmin()
    r2 = np.abs(lat-20).argmin() + 1
    r3 = np.abs(lon+70).argmin()
    r4 = np.abs(lon-40).argmin() + 1
    from scipy import stats
    import glob
    # point a to control directory, point b to forced directory, leaving '/*' at the end of the name
    a = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_577/atmos/item16222_monthly_mean/*')
    b = glob.glob('/network/aopp/hera/mad/bakerh/RPM/cpdn_extract_scripts-master/extracted_data/batch_578/atmos/item16222_monthly_mean/*')
    lst = np.arange(1, len(a)+1)
    anom = np.zeros((len(a), len(lat), len(lon)))
    NAOinds = np.zeros((2, len(a)))
    for i, item in enumerate(a):
        # load patches, point this to patch directory, leaving 'lst[i]' in the middle to loop over patch numbers
        anom[i, :] = ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert' + lst[i] + '_c040926.nc', 'SST')[0, :] - ncread('/home/bakerh/Documents/DPhil/CPDN/RPM/SSTfiles/sst_HadOIBl_bc_N96_clim_c040926_25months.nc', 'SST_cpl')[0, :]
        # load mslp from control run, alter str to variable name
        mslp_control = ncread(a[i], 'item16222_monthly_mean')[:, 0, :]/100
        # load mslp from forced run, alter str to variable name
        mslp = ncread(b[i], 'item16222_monthly_mean')[:, 0, :]/100
        # compute seasonal control MSLP fields
        mslp_cw, mslp_cs = win_sum(mslp_control)
        # compute forced seasonal MSLP fields
        mslp_w, mslp_s = win_sum(mslp)
        # compute MSLP projected onto EOF to give NAO index
        nao_cw = eof_response(np.expand_dims(mslp_cw[r1:r2, r3:r4], axis=0), eof[0, r1:r2, r3:r4], lat[r1:r2], lon[r3:r4])
        nao_cs = eof_response(np.expand_dims(mslp_cs[r1:r2, r3:r4], axis=0), eof[1, r1:r2, r3:r4], lat[r1:r2], lon[r3:r4])
        nao_w = eof_response(np.expand_dims(mslp_w[r1:r2, r3:r4], axis=0), eof[0, r1:r2, r3:r4], lat[r1:r2], lon[r3:r4])
        nao_s = eof_response(np.expand_dims(mslp_s[r1:r2, r3:r4], axis=0), eof[1, r1:r2, r3:r4], lat[r1:r2], lon[r3:r4])
        # compute EOF responses
        NAOinds[0, i] = nao_w - nao_cw
        NAOinds[1, i] = nao_s - nao_cs
        print('Done: ' + str(i+1) + ' out of ' + str(len(a)))
    # set SST anoms to -180 to 180, 90 to -90
    anom = np.concatenate((anom[:, ::-1, int(len(lon)/2):], anom[:, ::-1, :int(len(lon)/2)]), axis=2)
    regmap_smoothed = np.zeros((2, len(lat), len(lon)))
    regmap_smoothed_sig = np.zeros((2, len(lat), len(lon)))
    regmap_smoothed_sig[:] = np.NAN
    # perform linear regression
    for i in range(len(lat)):
        for j in range(len(lon)):
            for a in range(2):
                regmap_smoothed[a, i, j] = np.cov(anom[:, i, j], NAOinds[a, :])[0, 1]*3/(4*6371**2 * (lat[0]-lat[1]) * (lon[1]-lon[0]) * (np.pi/180)**2)
                r = np.cov(anom[:, i, j], NAOinds[a, :])[0, 1]/np.sqrt(np.cov(anom[:, i, j], NAOinds[a, :])[1, 1]*np.cov(anom[:, i, j], NAOinds[a, :])[0, 0])
                t = r*np.sqrt((len(a)-2)/(1-r**2))
                p = 1 - stats.norm.cdf(np.abs(t))
                sig = np.greater_equal(5, p*100*2).astype(int)
                if sig == 1:
                    regmap_smoothed_sig[a, i, j] = regmap_smoothed[a, i, j]
    return regmap_smoothed, regmap_smoothed_sig


def compute_reconstruction(gto, nao_20thc, lat, lon, region='glob'):
    # computes NAO reconstruction
    from scipy import interpolate
    from mpl_toolkits.basemap import shiftgrid
    lat_sst = np.linspace(89.5, -89.5, 180)
    lon_sst_shift = np.linspace(.5, 359.5, 360)
    # create seasonal indices
    sindices = np.zeros((147*3))
    windices = np.zeros((147*3))
    for i in range(147):
        sindices[3*i:3*(i+1)] = [17+12*i, 18+12*i, 19+12*i]
        windices[3*i:3*(i+1)] = [11+12*i, 12+12*i, 13+12*i]
    sindices = sindices.astype(int)
    windices = windices.astype(int)

    gto_interp = np.zeros((2, 180, 360))
    # point this at direcotry containing HadISST dataset
    sst = ncread('/network/aopp/hera/mad/bakerh/Reanalyses/HadISST/HadISST_sst.nc','sst')[:12*148]
    # remove ice
    sst_0 = np.where(sst < -999, 0, sst)
    # remove land
    gto_0 = np.where(np.isnan(gto), 0, gto)
    # put gto on same grid as SSTs
    for i in range(2):
        f = interpolate.interp2d(lon, lat[::-1], gto_0[i, :])
        gto_interp[i, :] = f(lon_sst_shift, lat_sst)
        gto_interp[i], lon1 = shiftgrid(181., gto_interp[i], lon_sst_shift,
                                        start=False)
    # specify region
    if region == 'atl':
        gto_interp[:, :, :90] = 0
        gto_interp[:, :, 210:] = 0
    elif region == 'ind':
        gto_interp[:, :, :210] = 0
        gto_interp[:, :, 290:] = 0
    elif region == 'pac':
        gto_interp[:, :, 90:290] = 0
    # weight SST grid and then project onto GTO
    meshlat = np.zeros([180, 360])
    meshlat[:, :] = np.expand_dims(lat_sst, axis=1)
    meshlatweight = np.cos(meshlat * np.pi/180) * 6371**2 * 1 * 1 * (np.pi/180)**2
    nao_recon_monthly = np.zeros((2, 3*147))
    nao_recon = np.zeros((2, 147))

    sst_w = sst_0[windices]
    sst_s = sst_0[sindices]
    sst_w = sst_w * meshlatweight
    sst_s = sst_s * meshlatweight
    # reconstruct NAO for each month
    for y in range(147*3):
        nao_recon_monthly[0, y] = np.nansum(sst_w[y]*gto_interp[0])
        nao_recon_monthly[1, y] = np.nansum((sst_s[y])*gto_interp[1])
    # mean over season
    for y in range(147):
        nao_recon[0, y] = np.mean(nao_recon_monthly[0, 3*y:3*(y+1)])
        nao_recon[1, y] = np.mean(nao_recon_monthly[1, 3*y:3*(y+1)])

    # now calibrate nao to observed nao
    calib_coef = calibration(nao_recon, nao_20thc)
    nao_recon = calibrate(nao_recon, calib_coef)
    return nao_recon


########################
# AUXILIARY FUNCTIONS  #
########################


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


def eofweight(unweighted, lat, lon):
    '''
    Outputs weighted data for projections (i.e. sqrt cos lat)

    Parameters
    ----------
    unweighted: array
        unweighted data
    lat: array
        latitudes
    lon: array
        lon values

    Outputs
    -------
    weighted: array
        weighted data
    '''
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.sqrt(np.cos(meshlat.transpose() * np.pi/180))
    weighted = unweighted * meshlatweight
    return weighted


def eof_response(response, eofn, lat, lon, corr='no'):
    '''
    Projects response onto eofn

    Parameters
    ----------
    response: array
        data to project on to eofn
    eofn: array
        nth eof
    lat: array
        latitude of data
    lon: array
        longitude of data

    Returns
    -------
    projection: array
        response projected on eofn
    '''
    eof_region = eofn.copy()
    r_region = response.copy()
    responsew = eofweight(r_region, lat, lon)
    eofn_w = eofweight(eof_region.copy(), lat, lon)
    projection = np.zeros(len(response))
    if corr == 'yes':
        for i in range(len(response)):
            projection[i] = (np.sum(responsew[i] * eofn_w) /
                             np.sqrt(np.sum(eofn_w*eofn_w) *
                                     np.sum(responsew[i]*responsew[i])))
    else:
        for i in range(len(response)):
            projection[i] = (np.sum(responsew[i] * eofn_w) /
                             (np.sum(eofn_w*eofn_w)))
    return projection


def eof_clim(clim, lat, lon, neof=1, region=[8, 57, 48, 118]):
    def eof_svd(clim, lat, lon, neof):
        '''
        EOFs of data

        Parameters
        ----------
        clim: array
            data to determine EOFs for. time x nlat x nlon

        Returns
        -------
        eofs: array
            First n EOFs. neof x nsigma x nlat
        var: array
            amount of variance explained by each EOF. neof x 3 array:
            var[:, 0] = the 1st neof eigenvalues of the covariance matrix
            var[:, 1] = the 1st neof explained variances (in %).
            var[:, 2] = the 1st neof errors (in %).
        pcs: array
            principal components. time x neof
        '''
        neof -= 1
        from numpy import linalg
        field = clim.copy()[:, region[0]:region[1], region[2]:region[3]]
        dim = field.shape
        field = eofweight(field, lat[region[0]:region[1]],
                          lon[region[2]:region[3]])
        aver = np.mean(field, axis=0)
        for i in range(dim[0]):
            field[i, :] = field[i, :] - aver
        field = np.reshape(field, (dim[0], dim[1] * dim[2]), order='F')

        if dim[0] > dim[1] * dim[2]:
            field = np.transpose(field)

        u, s, v = linalg.svd(field)
        v = np.transpose(v)
        s = s * s
        eigs = np.copy(s)
        s = s / np.sum(s)
        if dim[0] < dim[1] * dim[2]:
            u, v = v, u

        eofs = u[:, neof]
        pcs = v[:, neof]
        var = np.zeros(3)
        var[0] = eigs[neof]
        var[1] = s[neof] * 100
        var[2] = var[1] * np.sqrt(2/len(clim))
        eofs = np.transpose(eofs)
        eofs = np.reshape(eofs, (dim[1], dim[2]), order='F')
        return eofs, pcs, var

    def eof_regress(pcs, clim):
        '''
        Regresses original field on nth pcs time series

        Parameters
        ----------
        control: array
            original data to regress. time x nsigma x nlat
        eofn: int
            eof number to calculate
        pcs: array
            first principal component time series

        Returns
        -------
        eofn: array
            nth eof from field regressed onto pcsn
        '''
        from scipy import stats
        pcsn = pcs[:]
        pcsn = (pcsn - np.mean(pcsn)) / np.std(pcsn)
        eofn = np.zeros((np.ma.size(clim, axis=1),
                        np.ma.size(clim, axis=2)))
        for a in range(np.ma.size(clim, axis=1)):
            for b in range(np.ma.size(clim, axis=2)):
                eofn[a, b] = stats.linregress(pcsn, clim[:, a, b])[0]
        return eofn  # this is not normalized!

    # compute eofs
    eofs_un, pcs, var = eof_svd(clim, lat, lon, neof)
    # regress control onto pcs to get sensible units
    eofs = eof_regress(pcs, clim.copy())
    return eofs, var


def calibration(recon, ncar):
    # computes calibration coefficients for reconstructions
    from scipy import stats
    ncar = np.reshape(ncar, (-1, np.shape(ncar)[-1]))
    recon = np.reshape(recon, (-1, np.shape(recon)[-1]))
    nind = np.shape(recon)[0]
    conv = np.zeros((nind, 2))
    for i in range(nind):
        conv[i, :] = stats.linregress(recon[i, 1:-5], ncar[i])[:2]
    return conv


def calibrate(recon, conv):
    # calibrates values
    oshape = np.shape(recon)
    recon = np.reshape(recon, (-1, np.shape(recon)[-1]))
    nind = np.shape(recon)[0]
    recon_corrected = np.zeros((nind, np.shape(recon)[1]))
    for i in range(nind):
        recon_corrected[i, :] = recon[i]*conv[i, 0] + conv[i, 1]
    recon_corrected = np.reshape(recon_corrected, oshape)
    return recon_corrected
