#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 12:18:14 2018

@author: bakerh
"""

import numpy as np


def monthlymeans(item, code, level, region=[0, 144, 0, 192]):
    import glob
    from netCDF4 import Dataset
    sindices = np.zeros((10))
    for i in range(10):
        sindices[i] = i*12
    sindices = sindices.astype(int)
    output = {}

    exps = ['Plus15-Future_LCO2', 'Plus15-Future']
    exps = ['All-Nat', 'SST-Nat', 'GHG-Nat']

    for x, exp in enumerate(exps):
        a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/' +
                             exp + '/mon/' + item + '/*'))
        sdata = np.zeros((120, region[1]-region[0], region[3]-region[2]))
        tbar = np.zeros((12, region[1]-region[0], region[3]-region[2]))
        for e, d in enumerate(a):
            nc_fid = Dataset(d, 'r')
            sdata += nc_fid.variables[code][:120, level, region[0]:region[1],
                                              region[2]:region[3]]
            print('Done ' + str(item) + ' ' + str(exp) + ' ' + str(e+1))
        sdata /= len(a)
        for i in range(12):
            tbar[i] = np.mean(sdata[sindices+i], axis=(0))
        output[exp] = tbar#*86400
    return output


def w5rem(v_field):
    v_field_rem = np.zeros((12, 144, 192))
    for i in range(12):
        v_field_fft = np.fft.rfft2(v_field[i])
        v_field_fft[:, 2:] = 0
        v_field_rem[i] = np.fft.irfft2(v_field_fft)
    return v_field_rem


def rws_fprime(ua200, va200):
    from windspharm.standard import VectorWind
    from netcdfread import ncread, ncsave
    from scipy import interpolate

    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')   

    u = ua200['All-Nat']#[1]
    v = va200['All-Nat']#[1]
    u_sst = ua200['GHG-Nat']#[1]
    v_sst = va200['GHG-Nat']#[1]
    u_ghg = ua200['SST-Nat']#[1]
    v_ghg = va200['SST-Nat']#[1]

    v = w5rem(v)
    u = w5rem(u)

    uwnd = np.zeros((64, 128, len(u)))
    vwnd = np.zeros((64, 128, len(v)))
    uwnd_sst = np.zeros((64, 128, len(u_sst)))
    vwnd_sst = np.zeros((64, 128, len(u_sst)))
    uwnd_ghg = np.zeros((64, 128, len(u_ghg)))
    vwnd_ghg = np.zeros((64, 128, len(u_ghg)))

    u = np.transpose(u, (1, 2, 0))
    v = np.transpose(v, (1, 2, 0))
    u_sst = np.transpose(u_sst, (1, 2, 0))
    v_sst = np.transpose(v_sst, (1, 2, 0))
    u_ghg = np.transpose(u_ghg, (1, 2, 0))
    v_ghg = np.transpose(v_ghg, (1, 2, 0))

    for i in range(np.ma.size(uwnd, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u[:, :, i])
        uwnd[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v[:, :, i])
        vwnd[:, :, i] = g(lon42, lat42[::-1])
    for i in range(np.ma.size(uwnd_sst, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u_sst[:, :, i])
        uwnd_sst[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v_sst[:, :, i])
        vwnd_sst[:, :, i] = g(lon42, lat42[::-1])
    for i in range(np.ma.size(uwnd_ghg, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u_ghg[:, :, i])
        uwnd_ghg[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v_ghg[:, :, i])
        vwnd_ghg[:, :, i] = g(lon42, lat42[::-1])

    w = VectorWind(uwnd, vwnd)
    w_sst = VectorWind(uwnd_sst, vwnd_sst)
    w_ghg = VectorWind(uwnd_ghg, vwnd_ghg)

    eta = w.absolutevorticity()
    eta_sst = w_sst.absolutevorticity()
    eta_ghg = w_ghg.absolutevorticity()

    div = w.divergence()
    div_sst = w_sst.divergence()
    div_ghg = w_ghg.divergence()

    uchi, vchi = w.irrotationalcomponent()
    uchi_sst, vchi_sst = w_sst.irrotationalcomponent()
    uchi_ghg, vchi_ghg = w_ghg.irrotationalcomponent()

    etax, etay = w.gradient(eta)
    # etax_sst, etay_sst = w_sst.gradient(eta_sst)
    # etax_ghg, etay_ghg = w_ghg.gradient(eta_ghg)

    eta = np.transpose(eta, (2, 0, 1))
    eta_sst = np.transpose(eta_sst, (2, 0, 1))
    eta_ghg = np.transpose(eta_ghg, (2, 0, 1))

    div = np.transpose(div, (2, 0, 1))
    div_sst = np.transpose(div_sst, (2, 0, 1))
    div_ghg = np.transpose(div_ghg, (2, 0, 1))

    etax = np.transpose(etax, (2, 0, 1))
    etay = np.transpose(etay, (2, 0, 1))

    uchi = np.transpose(uchi, (2, 0, 1))
    uchi_sst = np.transpose(uchi_sst, (2, 0, 1))
    uchi_ghg = np.transpose(uchi_ghg, (2, 0, 1))
    vchi = np.transpose(vchi, (2, 0, 1))
    vchi_sst = np.transpose(vchi_sst, (2, 0, 1))
    vchi_ghg = np.transpose(vchi_ghg, (2, 0, 1))

    #f_ghg = -eta.mean(axis=0)*(div_ghg-div.mean(axis=0))-(uchi_ghg-uchi.mean(axis=0))*etax.mean(axis=0)-(vchi_ghg-vchi.mean(axis=0))*etay.mean(axis=0)
    #f_sst = -eta.mean(axis=0)*(div_sst-div.mean(axis=0))-(uchi_sst-uchi.mean(axis=0))*etax.mean(axis=0)-(vchi_sst-vchi.mean(axis=0))*etay.mean(axis=0)

    f_ghg = -eta*(div_ghg-div)-((uchi_ghg-uchi)*etax+(vchi_ghg-vchi)*etay)
    f_sst = -eta*(div_sst-div)-((uchi_sst-uchi)*etax+(vchi_sst-vchi)*etay)

    meshlon, meshlat = np.meshgrid(lon42, lat42)
    '''
    ncsave('/home/bakerh/Downloads/vort200_control', lat42, lon42, eta.mean(axis=0)-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/vort200_sst', lat42, lon42, eta_sst.mean(axis=0)-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/vort200_ghg', lat42, lon42, eta_ghg.mean(axis=0)-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/forcing_ghg', lat42, lon42, f_ghg.mean(axis=0), 'forcing')
    ncsave('/home/bakerh/Downloads/forcing_sst', lat42, lon42, f_sst.mean(axis=0), 'forcing')
    '''
    month = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ncsave('/home/bakerh/Downloads/vort200_control', month, lat42, lon42, eta-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/vort200_sst', month, lat42, lon42, eta_sst-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/vort200_ghg', month, lat42, lon42, eta_ghg-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/forcing_ghg', month, lat42, lon42, f_ghg, 'forcing')
    ncsave('/home/bakerh/Downloads/forcing_sst', month, lat42, lon42, f_sst, 'forcing')


def rws_fprime_co2(ua200, va200):
    from windspharm.standard import VectorWind
    from netcdfread import ncread, ncsave
    from scipy import interpolate

    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc', 'lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc', 'lat')
    # HadAM3P
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')   
    # MIROC5
    # lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/MIROC5/All-Hist/day/tas/tas_Aday_MIROC5_All-Hist_est1_v2-0_run001_20060101-20161231.nc','lat')
    # lon = np.arange(0, 360, 360/256)
    # CAM4
    #lat = np.linspace(-90, 90, 96)
    #lon = np.linspace(0, 357.5, 144)
    u = ua200['Plus15-Future_LCO2']
    v = va200['Plus15-Future_LCO2']
    u_ghg = ua200['Plus15-Future_HCO2']
    v_ghg = va200['Plus15-Future_HCO2']

    uwnd = np.zeros((64, 128, len(u)))
    vwnd = np.zeros((64, 128, len(v)))
    uwnd_ghg = np.zeros((64, 128, len(u_ghg)))
    vwnd_ghg = np.zeros((64, 128, len(u_ghg)))

    u = np.transpose(u, (1, 2, 0))
    v = np.transpose(v, (1, 2, 0))
    u_ghg = np.transpose(u_ghg, (1, 2, 0))
    v_ghg = np.transpose(v_ghg, (1, 2, 0))

    for i in range(np.ma.size(uwnd, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u[:, :, i])
        uwnd[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v[:, :, i])
        vwnd[:, :, i] = g(lon42, lat42[::-1])
    for i in range(np.ma.size(uwnd_ghg, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u_ghg[:, :, i])
        uwnd_ghg[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v_ghg[:, :, i])
        vwnd_ghg[:, :, i] = g(lon42, lat42[::-1])

    w = VectorWind(uwnd, vwnd)
    w_ghg = VectorWind(uwnd_ghg, vwnd_ghg)

    eta = w.absolutevorticity()
    eta_ghg = w_ghg.absolutevorticity()

    div = w.divergence()
    div_ghg = w_ghg.divergence()

    uchi, vchi = w.irrotationalcomponent()
    uchi_ghg, vchi_ghg = w_ghg.irrotationalcomponent()

    etax, etay = w.gradient(eta)
    # etax_ghg, etay_ghg = w_ghg.gradient(eta_ghg)

    eta = np.transpose(eta, (2, 0, 1))
    eta_ghg = np.transpose(eta_ghg, (2, 0, 1))

    div = np.transpose(div, (2, 0, 1))
    div_ghg = np.transpose(div_ghg, (2, 0, 1))

    etax = np.transpose(etax, (2, 0, 1))
    etay = np.transpose(etay, (2, 0, 1))

    uchi = np.transpose(uchi, (2, 0, 1))
    uchi_ghg = np.transpose(uchi_ghg, (2, 0, 1))
    vchi = np.transpose(vchi, (2, 0, 1))
    vchi_ghg = np.transpose(vchi_ghg, (2, 0, 1))

    f_ghg = -eta*(div_ghg-div)-((uchi_ghg-uchi)*etax+(vchi_ghg-vchi)*etay)

    meshlon, meshlat = np.meshgrid(lon42, lat42)
    '''
    ncsave('/home/bakerh/Downloads/vort200_control', lat42, lon42, eta.mean(axis=0)-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/forcing_ghg', lat42, lon42, f_ghg.mean(axis=0), 'forcing')
    '''
    month = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
    ncsave('/home/bakerh/Downloads/vort200_control', month, lat42, lon42, eta-2*np.sin(meshlat*np.pi/180)*7.2921e-5, 'vorticity')
    ncsave('/home/bakerh/Downloads/forcing_ghg', month, lat42, lon42, f_ghg, 'forcing')


def windprops(ua, va, lat, lon, exps=['All-Nat', 'SST-Nat']):
    from windspharm.standard import VectorWind

    def diverg(u, v, lat, lon):
        meshlon, meshlat = np.meshgrid(lon, lat)
        phi = meshlat*np.pi/180
        cs_phi = np.cos(phi)
        dphi = lat*np.pi/180
        dtheta = lon*np.pi/180
        a = 6.371e6

        v_y = np.gradient(v*cs_phi, dphi, axis=1) / (a*cs_phi)
        u_x = np.gradient(u, dtheta, axis=2) / (a*cs_phi)
        div = u_x + v_y
        return div

    u = ua[exps[0]]
    v = va[exps[0]]
    u_ghg = ua[exps[1]]
    v_ghg = va[exps[1]]

    div = diverg(u, v, lat, lon)
    div_ghg = diverg(u_ghg, v_ghg, lat, lon)

    u = np.transpose(u, (1, 2, 0))
    v = np.transpose(v, (1, 2, 0))
    u_ghg = np.transpose(u_ghg, (1, 2, 0))
    v_ghg = np.transpose(v_ghg, (1, 2, 0))

    w = VectorWind(u, v)
    w_ghg = VectorWind(u_ghg, v_ghg)
    vp = w.velocitypotential()
    vp_ghg = w_ghg.velocitypotential()

    vp = np.transpose(vp, (2, 0, 1))
    vp_ghg = np.transpose(vp_ghg, (2, 0, 1))

    delta_vp = vp_ghg - vp
    delta_div = div_ghg - div

    return div, delta_div, vp, delta_vp


def ming_forcing(ua, va, psl):
    from windspharm.standard import VectorWind
    from netcdfread import ncread, ncsave
    from scipy import interpolate

    lon42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat42 = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/Plus15-Future_LCO2/day/ua/item15201_daily_mean_a00b_2090-01_2100-12.nc', 'latitude1')   
    lat145 = ncread('/network/aopp/hera/mad/bakerh/HAPPI/HadAM3P-N96/All-Hist/mon/tas/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')

    u = ua['All-Nat'][1]
    v = va['All-Nat'][1]
    u_sst = ua['GHG-Nat'][1]
    v_sst = va['GHG-Nat'][1]
    u_ghg = ua['SST-Nat'][1]
    v_ghg = va['SST-Nat'][1]
    p = psl['All-Nat'][1]
    p_sst = psl['GHG-Nat'][1]
    p_ghg = psl['SST-Nat'][1]

    uwnd = np.zeros((64, 128, len(u)))
    vwnd = np.zeros((64, 128, len(v)))
    uwnd_sst = np.zeros((64, 128, len(u_sst)))
    vwnd_sst = np.zeros((64, 128, len(u_sst)))
    uwnd_ghg = np.zeros((64, 128, len(u_ghg)))
    vwnd_ghg = np.zeros((64, 128, len(u_ghg)))
    ps = np.zeros((64, 128, len(p)))
    ps_sst = np.zeros((64, 128, len(p_sst)))
    ps_ghg = np.zeros((64, 128, len(p_ghg)))

    u = np.transpose(u, (1, 2, 0))
    v = np.transpose(v, (1, 2, 0))
    u_sst = np.transpose(u_sst, (1, 2, 0))
    v_sst = np.transpose(v_sst, (1, 2, 0))
    u_ghg = np.transpose(u_ghg, (1, 2, 0))
    v_ghg = np.transpose(v_ghg, (1, 2, 0))
    p = np.transpose(p, (1, 2, 0))
    p_sst = np.transpose(p_sst, (1, 2, 0))
    p_ghg = np.transpose(p_ghg, (1, 2, 0))

    for i in range(np.ma.size(uwnd, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u[:, :, i])
        uwnd[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v[:, :, i])
        vwnd[:, :, i] = g(lon42, lat42[::-1])
        j = interpolate.interp2d(lon, lat145[::-1], p[:, :, i])
        ps[:, :, i] = j(lon42, lat42[::-1])
    for i in range(np.ma.size(uwnd_sst, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u_sst[:, :, i])
        uwnd_sst[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v_sst[:, :, i])
        vwnd_sst[:, :, i] = g(lon42, lat42[::-1])
        j = interpolate.interp2d(lon, lat145[::-1], p_sst[:, :, i])
        ps_sst[:, :, i] = j(lon42, lat42[::-1])
    for i in range(np.ma.size(uwnd_ghg, axis=2)):
        h = interpolate.interp2d(lon, lat[::-1], u_ghg[:, :, i])
        uwnd_ghg[:, :, i] = h(lon42, lat42[::-1])
        g = interpolate.interp2d(lon, lat[::-1], v_ghg[:, :, i])
        vwnd_ghg[:, :, i] = g(lon42, lat42[::-1])
        j = interpolate.interp2d(lon, lat145[::-1], p_ghg[:, :, i])
        ps_ghg[:, :, i] = j(lon42, lat42[::-1])

    w = VectorWind(uwnd, vwnd)
    w_sst = VectorWind(uwnd_sst, vwnd_sst)
    w_ghg = VectorWind(uwnd_ghg, vwnd_ghg)

    psx, psy = w.gradient(ps)
    psx_sst, psy_sst = w_sst.gradient(ps_sst)
    psx_ghg, psy_ghg = w_ghg.gradient(ps_ghg)

    omega = np.mean(uwnd * psx + vwnd * psy, axis=2)
    omega_sst = uwnd_sst * psx_sst + vwnd_sst * psy_sst
    omega_ghg = uwnd_ghg * psx_ghg + vwnd_ghg * psy_ghg

    omega_sst = np.transpose(omega_sst, (2, 0, 1))
    omega_ghg = np.transpose(omega_ghg, (2, 0, 1))

    oforc_sst = np.mean(omega_sst-omega, axis=0)
    oforc_ghg = np.mean(omega_ghg-omega, axis=0)

    return omega, oforc_sst, oforc_ghg


def sens():
    import numpy as np
    from netcdfread import ncread
    import glob
    a = sorted(glob.glob('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/sens/*'))
    v = np.zeros((36, 26, 64, 128))
    for i in range(36):
        v[i] = ncread(a[i], 'V')
    return v


def response(v_model, model='ghg', month='J', r=[0, 90, 0, 360], ln=14):
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

    def eof_response(response, eofn, lati, long, region=[-90, 90, 0, 360],
                     corr='yes'):
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
        lat = lati[(region[0] <= lati) & (lati <= region[1])]
        lon = long[(region[2] <= long) & (long <= region[3])]
        eof_region = eofn.copy()[(region[0] <= lati) & (lati <= region[1]), :]
        eof_region = eof_region[:, (region[2] <= long) & (long <= region[3])]
        r_region = response.copy()[:, (region[0] <= lati) & (lati <= region[1]), :]
        r_region = r_region[:, :, (region[2] <= long) & (long <= region[3])]
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
                                 np.sqrt(np.sum(eofn_w*eofn_w)))
        return projection

    import numpy as np
    from netcdfread import ncread
    import glob
    lon = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lati = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')
    lat = np.copy(lati[::-1])
    a = glob.glob('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/outputs/' +
                  model + '/' + month + '/*')
    v_f = np.zeros((32*ln, 21, 64, 128))
    for i in range(32*ln):
        v_f[i] = ncread(a[i], 'V') * 1e6
    v_f = np.mean(v_f[:, 16:], axis=1) - v_f[:, 0]
    w = eof_response(v_f, v_model, lati, lon, region=r,
                     corr='no')
    f_per = eof_response(v_f, v_model, lati, lon, region=r)
    w = np.reshape(w, (ln, 32))
    f_p = np.reshape(f_per, (ln, 32))
    weights = np.zeros((64, 128))
    forcing_p = np.zeros((64, 128))

    for LATraster in range(1, ln+1):
        for LONraster in range(1, 33):
            w1 = w[LATraster-1, LONraster-1]
            f1 = f_p[LATraster-1, LONraster-1]
            if np.mod(LATraster, 2) == 0:
                LONraster = LONraster + .5
            if np.mod(LATraster, 2) != 0 and LONraster == 1:
                for ii in range(int(LONraster*4-4), int(LONraster*4-4+2)):
                    for jj in range(LATraster+32-2, LATraster+32+3):
                        weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*4-4)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
                        forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*4-4)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
                for ii in range(126, 128):
                    for jj in range(LATraster+32-2, LATraster+32+3):
                        weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-360)/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
                        forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-360)/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
            elif np.mod(LATraster, 2) == 0 and LONraster == 32.5:
                for ii in range(int(LONraster*4-4-2), int(LONraster*4-4+2)):
                    for jj in range(LATraster+32-2, LATraster+32+3):
                        weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*4-4)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
                        forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*4-4)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
            else:
                for ii in range(int(LONraster*4-4-2), int(LONraster*4-4+3)):
                    for jj in range(LATraster+32-2, LATraster+32+3):
                        weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*4-4)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)
                        forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*4-4)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[LATraster+32])/(lat[31]-lat[33])))**2)

    weights = weights[::-1, :]
    forcing_p = forcing_p[::-1, :] * 2.5e-11
    v_perfect = np.mean(v_f.transpose(1, 2, 0)*f_per, axis=2)
    return weights, forcing_p, v_perfect, w, f_per


def response_1(v_model, model='ghg', month='J', r=[0, 90, 0, 360]):
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

    def eof_response(response, eofn, lati, long, region=[-90, 90, 0, 360],
                     weigh='yes'):
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
        lat = lati[(region[0] <= lati) & (lati <= region[1])]
        lon = long[(region[2] <= long) & (long <= region[3])]
        eof_region = eofn.copy()[(region[0] <= lati) & (lati <= region[1]), :]
        eof_region = eof_region[:, (region[2] <= long) & (long <= region[3])]
        r_region = response.copy()[:, (region[0] <= lati) & (lati <= region[1]), :]
        r_region = r_region[:, :, (region[2] <= long) & (long <= region[3])]
        responsew = eofweight(r_region, lat, lon)
        eofn_w = eofweight(eof_region.copy(), lat, lon)
        projection = np.zeros(len(response))
        if weigh == 'yes':
            for i in range(len(response)):
                projection[i] = (np.sum(responsew[i] * eofn_w) /
                                 (np.sum(eofn_w*eofn_w)))
        else:
            for i in range(len(response)):
                projection[i] = (np.sum(responsew[i] * eofn_w) /
                                 ((np.sum(responsew[i]*responsew[i]))))
        return projection

    import numpy as np
    from netcdfread import ncread
    import glob
    lon = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lati = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')
    lat = np.copy(lati[::-1])
    a = glob.glob('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/outputs/single/' +
                  model + '/' + month + '/*')
    v_f = np.zeros((64, 21, 64, 128))
    for i in range(64):
        v_f[i] = ncread(a[i], 'V')
    v_f = (np.mean(v_f[:, 16:], axis=1) - v_f[:, 0]) * 1e6
    w = eof_response(v_f, v_model, lati, lon, region=r,
                     weigh='yes')
    f_per = eof_response(v_f, v_model, lati, lon, region=r, weigh='no')
    weights = np.zeros((64, 128))
    forcing_p = np.zeros((64, 128))

    for LONraster in range(1, 65):
        w1 = w[LONraster-1]
        f1 = f_per[LONraster-1]
        if LONraster == 1:
            for ii in range(int(LONraster*2-2), int(LONraster*2-2+3)):
                for jj in range(34-2, 34+3):
                    weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*2-2)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
                    forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*2-2)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
            for ii in range(126, 128):
                for jj in range(34-2, 34+3):
                    weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-360)/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
                    forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-360)/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
        elif LONraster == 64:
            for ii in range(int(LONraster*2-2-2), int(LONraster*2-2+2)):
                for jj in range(34-2, 34+3):
                    weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*2-2)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
                    forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*2-2)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
        else:
            for ii in range(int(LONraster*2-2-2), int(LONraster*2-2+3)):
                for jj in range(34-2, 34+3):
                    weights[jj, ii] += w1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*2-2)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)
                    forcing_p[jj, ii] += f1 * ((np.cos(0.5*np.pi*(lon[int(ii)]-lon[int(LONraster*2-2)])/(lon[3]-lon[1])))**2)*((np.cos(0.5*np.pi*(lat[jj]-lat[34])/(lat[31]-lat[33])))**2)

    weights = weights[::-1, :]
    forcing_p = forcing_p[::-1, :] * 2.5e-11
    v_perfect = np.sum(v_f.transpose(1, 2, 0)*f_per, axis=2)
    return weights, forcing_p, v_perfect, w, f_per


def maplotter(pdata, colormax=1, colormin=-999, title='', output='yes'):
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
    lon = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')
    pdata, lon = shiftgrid(180., pdata, lon, start=False)
    pdata, lon = addcyclic(pdata, lon)
    plt.figure()
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
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(colormin, colormax, 17)
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    b = plt.colorbar(plot, orientation='horizontal',
                     aspect=50, shrink=0.75, spacing='proportional')
    b.set_label(label=r'V200 (1x10$^{-6}$ ms$^{-1}$)')
    #b.set_label(label=r'Forcing (1x10$^{-11}$ s$^{-2}$)')
    #b.set_label(label=r'Weights')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def maplot_ll(pdata, lat, lon, colormax=1, colormin=-999, title='', shift='yes', precip='no'):
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
    if colormin == -999:
        colormin = -colormax
    plt.figure()
    if shift == 'yes':
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
    b.set_label(label=r'EOF (K s.d.$^{-1}$)')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def maplot_llc(pdata, pdatac, lat, lon, colormax=1, colormin=-999, title='', precip='no'):
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
    if colormin == -999:
        colormin = -colormax
    plt.figure()
    pdata, lon1 = shiftgrid(180., pdata, lon, start=False)
    pdata, lon1 = addcyclic(pdata, lon1)

    pdatac, lon = shiftgrid(180., pdatac, lon, start=False)
    pdatac, lon = addcyclic(pdatac, lon)

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
    #ctrsc = np.concatenate((np.linspace(-40, -5, 8), np.linspace(5, 40, 8)))
    ctrsc = np.concatenate((np.linspace(-2, -0.25, 8), np.linspace(.25, 2, 8)))
    ctrsc = [0.05]
    plot = m.contourf(x, y, pdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    m.contour(x, y, pdatac, ctrsc, colors='k')
    b = plt.colorbar(plot, orientation='horizontal',
                     aspect=50, shrink=0.75, spacing='proportional')
    #b.set_label(label=r'u200 difference (ms$^{-1}$)')
    b.set_label(label='Divergence ($1x10^{-5}$ s$^{-1}$)')
    #b.set_label(label='Velocity potential ($1x10^{7}$ m$^2$s$^{-1}$)') 
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def animate_rw(data, colorlimit=0.5, ts=6, title='', save='no',
               location='/home/bakerh/Documents/DPhil/Figures/Python/BTmodel/animations/new.mp4'):
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
    lon = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lon')
    lat = ncread('/network/aopp/hera/mad/bakerh/BTmodel_COR/main/inputs/forcing_ghg.nc','lat')

    summer = np.copy(data[::ts])

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    summer, lon = shiftgrid(180., summer, lon, start=False)
    summer, lon = addcyclic(summer, lon)
    meshlon, meshlat = np.meshgrid(lon, lat)
    ctrs = np.linspace(-colorlimit, colorlimit, 17)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    cmp = newcmap

    def updatefig(i):
        fig.clear()

        ax1 = fig.add_subplot(1, 1, 1)
        m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                    llcrnrlon=-180, urcrnrlon=180, resolution='c')
        m.drawcoastlines()
        m.drawcountries()
        m.drawmapboundary(linewidth=2)
        x, y = m(meshlon, meshlat)
        plot = m.contourf(x, y, summer[i]-summer[0], ctrs,
                          cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                          extend='both')
        # parallels = m.drawparallels(np.arange(-90., 91., 15.))
        # meridians = m.drawmeridians(np.arange(-180., 181., 30.))
        # m.drawparallels(parallels, labels=[True, True, True, True])
        # m.drawmeridians(meridians, labels=[True, True, True, True])
        b = fig.colorbar(plot, aspect=50, shrink=0.75,
                         spacing='proportional', orientation='horizontal',
                         extend='max', pad=0.05)
        b.set_label(label=r'V200 (1x10$^{-6}$ ms$^{-1}$)')
        j = i/(24/ts)
        ax1.set_title(' t=' + str('%02.3f' % j), loc='left', y=.9)
        plt.subplots_adjust(top=.99, bottom=0.01,
                            left=.05, right=.95)
        plt.suptitle(title, y=.85)
        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig, len(summer))
    if save == 'yes':
        anim.save(location, codec='mpeg4', bitrate=8000, dpi=300)
    return anim


def project(fielda, fieldb, lat1, lon1, lat2, lon2):
    from scipy import interpolate
    lat = np.linspace(90, -90, 145)
    lon = np.arange(0, 360, 360/192)
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.sqrt(np.cos(meshlat.transpose() * np.pi/180))
    h = interpolate.interp2d(lon1, lat1[::-1], fielda)
    field1 = h(lon, lat[::-1])
    h = interpolate.interp2d(lon2, lat2[::-1], fieldb)
    field2 = h(lon, lat[::-1])
    field1 = (field1*meshlatweight)[16:49]
    field2 = (field2*meshlatweight)[16:49]
    #field1 = np.hstack((((field1 - (field1*meshlatweight**2)[60:77,:32].mean()/(meshlatweight[60:77,:32]**2).mean())*meshlatweight)[60:77,:32], ((field1 - (field1*meshlatweight**2)[60:77,160:].mean()/(meshlatweight[60:77,160:]**2).mean())*meshlatweight)[60:77,160:]))
    #field2 = np.hstack((((field2 - (field2*meshlatweight**2)[60:77,:32].mean()/(meshlatweight[60:77,:32]**2).mean())*meshlatweight)[60:77,:32], ((field2 - (field2*meshlatweight**2)[60:77,160:].mean()/(meshlatweight[60:77,160:]**2).mean())*meshlatweight)[60:77,160:]))
    p_corr = (np.sum(field1 * field2) / np.sqrt(np.sum(field1*field1) *
                                                np.sum(field2*field2)))
    return p_corr











































