# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 10:20:45 2017

@author: bakerh
"""
import numpy as np


def main(var, global_tas):
    from scipy import interpolate
    from netcdfread import ncread
    import glob
    sindices = np.zeros((95*3))
    for i in range(95):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    sindices = sindices.astype(int)
    lat_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'latitude')
    lon_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'longitude')

    concs = {}
    concs['inmcm4'] = 1084.44
    concs['MIROC-ESM'] = 533.57
    concs['CanESM2'] = 532.45
    concs['MPI-ESM-LR'] = 555.37
    concs['GFDL-ESM2M'] = 778.35
    concs['HadGEM2-CC'] = 759.01
    concs['bcc-csm1-1'] = 675.27
    concs['IPSL-CM5A-MR'] = 656.68
    concs['GFDL-ESM2G'] = 770.60
    concs['IPSL-CM5B-LR'] = 666.95
    concs['MIROC-ESM-CHEM'] = 537.12
    concs['IPSL-CM5A-LR'] = 604.38
    concs['HadGEM2-ES'] = 740.65

    yrs = {}
    # yrs['CESM1-BGC'] = 2019.6038556835267
    yrs['CanESM2'] = 2015.6501388153354
    yrs['GFDL-ESM2G'] = 2040.2726255644302
    yrs['GFDL-ESM2M'] = 2038.4461621682028
    yrs['HadGEM2-CC'] = 2030.5237032804
    yrs['HadGEM2-ES'] = 2028.1144760733512
    yrs['IPSL-CM5A-LR'] = 2013.7490131087695
    yrs['IPSL-CM5A-MR'] = 2018.568003728501
    yrs['IPSL-CM5B-LR'] = 2025.1511344273508
    yrs['MIROC-ESM'] = 2021.3187712711444
    yrs['MIROC-ESM-CHEM'] = 2018.89262095655
    yrs['MPI-ESM-LR'] = 2020.0421238126532
    # yrs['NorESM1-ME'] = 2034.3684277126238
    yrs['bcc-csm1-1'] = 2022.8345355891056
    yrs['inmcm4'] = 2045.194355545509

    yrs2 = {}

    # yrs2['CESM1-BGC'] = 2034.7262797142157,
    yrs2['CanESM2'] = 2028.6584559023813,
    yrs2['GFDL-ESM2G'] = 2056.5526316368,
    yrs2['GFDL-ESM2M'] = 2053.6514333101513,
    yrs2['HadGEM2-CC'] = 2043.093047747961,
    yrs2['HadGEM2-ES'] = 2039.480727397444,
    yrs2['IPSL-CM5A-LR'] = 2029.344424079118,
    yrs2['IPSL-CM5A-MR'] = 2032.351689422151,
    yrs2['IPSL-CM5B-LR'] = 2039.2780954952768,
    yrs2['MIROC-ESM'] = 2031.144339409215,
    yrs2['MIROC-ESM-CHEM'] = 2030.7223697396423,
    yrs2['MPI-ESM-LR'] = 2039.3393134034848,
    # yrs2['NorESM1-ME'] = 2047.8103280489413,
    yrs2['bcc-csm1-1'] = 2038.2589979997997,
    yrs2['inmcm4'] = 2058.6460082149406

    tcr = {}
    # tcr['CESM1-BGC'] = 1.7
    tcr['CanESM2'] = 2.4
    tcr['GFDL-ESM2G'] = 1.1
    tcr['GFDL-ESM2M'] = 1.3
    tcr['HadGEM2-CC'] = np.nan
    tcr['HadGEM2-ES'] = 2.5
    tcr['IPSL-CM5A-LR'] = 2
    tcr['IPSL-CM5A-MR'] = 2
    tcr['IPSL-CM5B-LR'] = 1.5
    tcr['MIROC-ESM'] = 2.2
    tcr['MIROC-ESM-CHEM'] = np.nan
    tcr['MPI-ESM-LR'] = 2.0
    # tcr['NorESM1-ME'] = 1.6
    tcr['bcc-csm1-1'] = 1.7
    tcr['inmcm4'] = 1.3

    tcre = {}
    # tcre['CESM1-BGC'] = 1.7
    tcre['CanESM2'] = 2.375
    tcre['GFDL-ESM2G'] = 0.825
    tcre['GFDL-ESM2M'] = 1.075
    tcre['HadGEM2-CC'] = np.nan
    tcre['HadGEM2-ES'] = 2.1
    tcre['IPSL-CM5A-LR'] = 1.6
    tcre['IPSL-CM5A-MR'] = 1.575
    tcre['IPSL-CM5B-LR'] = 1.2
    tcre['MIROC-ESM'] = 2.15
    tcre['MIROC-ESM-CHEM'] = np.nan
    tcre['MPI-ESM-LR'] = 1.6
    # tcre['NorESM1-ME'] = 1.6
    tcre['bcc-csm1-1'] = 1.4
    tcre['inmcm4'] = 1

    tx90p_n96 = {}
    tx90p_n96_2 = {}
    items = []
    tx90p_MID = []
    tx90p_TROP = []
    ts_tx90p_MID = []
    ts_tx90p_TROP = []
    emissions = []
    years = []
    tcrs = []
    tcres = []
    global_tass = []
    a = glob.glob('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/*')
    #a1 = a[0]
    #a[0] = a[1]
    #a[1] = a1
    for i in a:
        items.append(i[78+2*len(var)-10:-30])
    for i in items:
        print(i)
        if i == 'HadGEM2-CC':
            tx90p = np.array(ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-210012.nc', var + 'ETCCDI')[12:, :])
            lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-210012.nc', 'lat')
            lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-210012.nc', 'lon')
        elif i == 'HadGEM2-ES':
            tx90p = np.array(ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-229912.nc', var + 'ETCCDI')[12:1152, :])
            lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-229912.nc', 'lat')
            lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-229912.nc', 'lon')        
        elif i == 'bcc-csm1-1':
            tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
            lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lat')
            lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lon')
        elif i == 'IPSL-CM5A-LR':
            tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
            lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lat')
            lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lon')
        elif i == 'MPI-ESM-LR':
            tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
            lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lat')
            lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lon')
        else:
            tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-210012.nc', var + 'ETCCDI')[:1140, :]
            lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-210012.nc', 'lat')
            lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-210012.nc', 'lon')
        tx90p = tx90p[sindices]
        tx90p2 = np.copy(tx90p)
        ts_tx90p = np.copy(tx90p)
        # idx = (np.abs(conversion_0[:, 2]-concs[i])).argmin()
        # yr = conversion_0[idx, 0]
        yr = int(np.round(yrs[i]))
        yr2 = int(np.round(yrs2[i]))
        tx90p = np.mean(tx90p[3*(yr-2010):3*(yr-2000), :], axis=0)
        tx90p2 = np.mean(tx90p2[3*(yr2-2010):3*(yr2-2000), :], axis=0)
        f = interpolate.interp2d(lon, lat[::-1], tx90p)
        f2 = interpolate.interp2d(lon, lat[::-1], tx90p2)
        ts_tx90p = np.reshape(ts_tx90p, (95, 3, len(lat), len(lon)))
        ts_tx90p = np.mean(ts_tx90p, axis=1)
        ts_tx90p_mid = np.zeros((95))
        ts_tx90p_trop = np.zeros((95))
        '''
        for j in range(95):
            g = interpolate.interp2d(lon, lat[::-1], ts_tx90p[j])
            ts_tx90p_mid[j] = latweightmean(g(lon_n96, lat_n96), msk='MID')
            ts_tx90p_trop[j] = latweightmean(g(lon_n96, lat_n96), msk='TROP')
        ts_tx90p_MID.append(ts_tx90p_mid)
        ts_tx90p_TROP.append(ts_tx90p_trop)
        '''
        tx90p_n96[i] = f(lon_n96, lat_n96)
        tx90p_n96_2[i] = f2(lon_n96, lat_n96)
        emissions.append(concs[i])
        tcres.append(tcre[i])
        tcrs.append(tcr[i])
        years.append(yrs[i])
        global_tass.append(global_tas[i + '_r1i1p1'])
        tx90p_MID.append(latweightmean(f(lon_n96, lat_n96), msk='MID'))
        tx90p_TROP.append(latweightmean(f(lon_n96, lat_n96), msk='TROP'))
    if var == 'tx90p':
        d = 0.9
    else:
        d = 1
    t_l = (tx90p_n96['CanESM2'] + tx90p_n96['MIROC-ESM-CHEM'] +
           tx90p_n96['IPSL-CM5A-LR'])*d/3
    t_h = (tx90p_n96['GFDL-ESM2G'] + tx90p_n96['inmcm4'] +
           tx90p_n96['HadGEM2-CC'])*d/3
    t2 = (tx90p_n96_2['CanESM2'] + tx90p_n96_2['MIROC-ESM-CHEM'] +
           tx90p_n96_2['IPSL-CM5A-LR'] + tx90p_n96_2['GFDL-ESM2G'] +
           tx90p_n96_2['inmcm4'] + tx90p_n96_2['HadGEM2-CC'])*d/6
    t15 = (tx90p_n96['CanESM2'] + tx90p_n96['MIROC-ESM-CHEM'] +
           tx90p_n96['IPSL-CM5A-LR'] + tx90p_n96['GFDL-ESM2G'] +
           tx90p_n96['inmcm4'] + tx90p_n96['HadGEM2-CC'])*d/6
    '''
    for i in range(len(years)):
        ts_tx90p_MID[i] *= d
        ts_tx90p_TROP[i] *= d
        tx90p_MID[i] *= d
        tx90p_TROP[i] *= d
    '''
    return t_l, t_h, t15, t2, tcrs, tcres, emissions, years, global_tass#, ts_tx90p_MID, ts_tx90p_TROP


def comp(model, var='tx90p'):
    from scipy import interpolate
    from netcdfread import ncread
    sindices = np.zeros((95*3))
    for i in range(95):
        sindices[3*i:3*(i+1)] = [5+12*i, 6+12*i, 7+12*i]
    sindices = sindices.astype(int)
    lat_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'latitude')
    lon_n96 = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'longitude')

    yrs = {}
    # yrs['CESM1-BGC'] = 2019.6038556835267
    yrs['CanESM2'] = 2015.6501388153354
    yrs['GFDL-ESM2G'] = 2040.2726255644302
    yrs['GFDL-ESM2M'] = 2038.4461621682028
    yrs['HadGEM2-CC'] = 2030.5237032804
    yrs['HadGEM2-ES'] = 2028.1144760733512
    yrs['IPSL-CM5A-LR'] = 2013.7490131087695
    yrs['IPSL-CM5A-MR'] = 2018.568003728501
    yrs['IPSL-CM5B-LR'] = 2025.1511344273508
    yrs['MIROC-ESM'] = 2021.3187712711444
    yrs['MIROC-ESM-CHEM'] = 2018.89262095655
    yrs['MPI-ESM-LR'] = 2020.0421238126532
    # yrs['NorESM1-ME'] = 2034.3684277126238
    yrs['bcc-csm1-1'] = 2022.8345355891056
    yrs['inmcm4'] = 2045.194355545509

    yrs26 = {}
    yrs26['CanESM2'] = 2015.9168461250079
    yrs26['HadGEM2-ES'] = 2031.2511570532158
    yrs26['IPSL-CM5A-LR'] = 2013.6504936554309
    yrs26['IPSL-CM5A-MR'] = 2019.6423318791442
    yrs26['MIROC-ESM'] = 2022.8139042931136
    yrs26['MIROC-ESM-CHEM'] = 2016.8634829024122
    yrs26['MPI-ESM-LR'] = 2025.4963052174598
    # yrs['NorESM1-ME'] = 2044.9730257768917
    yrs26['bcc-csm1-1'] = 2025.1860768601402

    i = model
    if i == 'HadGEM2-CC':
        tx90p = np.array(ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-210012.nc', var + 'ETCCDI')[12:, :])
        tx90p26 = np.array(ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp26/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp26_r1i1p1_200501-210012.nc', var + 'ETCCDI')[12:, :])
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-210012.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-210012.nc', 'lon')
    elif i == 'HadGEM2-ES':
        tx90p = np.array(ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-229912.nc', var + 'ETCCDI')[12:1152, :])
        tx90p26 = np.array(ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp26/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp26_r1i1p1_200501-229912.nc', var + 'ETCCDI')[12:1152, :])
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-229912.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200501-229912.nc', 'lon')        
    elif i == 'bcc-csm1-1':
        tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
        tx90p26 = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp26/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp26_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]       
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lon')
    elif i == 'IPSL-CM5A-LR':
        tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
        tx90p26 = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp26/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp26_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lon')
    elif i == 'MPI-ESM-LR':
        tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
        tx90p26 = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp26/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp26_r1i1p1_200601-230012.nc', var + 'ETCCDI')[:1140, :]
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-230012.nc', 'lon')
    else:
        tx90p = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-210012.nc', var + 'ETCCDI')[:1140, :]
        tx90p26 = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp26/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp26_r1i1p1_200601-210012.nc', var + 'ETCCDI')[:1140, :]
        lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-210012.nc', 'lat')
        lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/climdex/CMIP5/rcp85/' + var + '/' + var + 'ETCCDI_mon_' + str(i) + '_rcp85_r1i1p1_200601-210012.nc', 'lon')
    tx90p = tx90p[sindices]
    tx90p26 = tx90p26[sindices]
    # idx = (np.abs(conversion_0[:, 2]-concs[i])).argmin()
    # yr = conversion_0[idx, 0]
    yr = np.round(yrs[i])
    yr26 = np.round(yrs26[i])
    tx90p = np.mean(tx90p[3*(yr-2010):3*(yr-2000), :], axis=0)
    tx90p26 = np.mean(tx90p26[3*(yr26-2010):3*(yr26-2000), :], axis=0)
    f = interpolate.interp2d(lon, lat[::-1], tx90p)
    f26 = interpolate.interp2d(lon, lat[::-1], tx90p26)
    tx90p_n96 = f(lon_n96, lat_n96)*.9
    tx90p26_n96 = f26(lon_n96, lat_n96)*.9
    plotall(tx90p26_n96, tx90p_n96, 20, pltlbl=model + ' RCP8.5 - RCP2.6 TX90p')


def plotall(data_l, data_h, colorlimit, mask='yes', precip='no', pltlbl='RCP8.5 TX90p'):
    from netcdfread import ncread
    import matplotlib as mpl
    from matplotlib import pyplot as plt
    from mpl_toolkits.basemap import Basemap, shiftgrid, addcyclic
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    summer = data_h - data_l
    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    if mask == 'yes':
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
        summer = np.ma.masked_array(summer, mask=np.logical_not(lsm))

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

    ax1 = fig.add_subplot(1, 1, 1)
    m = Basemap(projection='moll', lon_0=0, resolution='c')
    m.drawcoastlines()
    m.drawcountries()
    m.drawmapboundary(linewidth=2)
    x, y = m(meshlon, meshlat)
    plot = m.contourf(x, y, summer, ctrs,
                      cmap=cmp, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')

    hfont = {'fontname': 'Arial'}
    c = fig.colorbar(plot, ax=ax1, orientation='horizontal',
                                 spacing='proportional', aspect=50)
    c.set_label(label='Difference (days season$\mathregular{^{-1}}$)', size=20)
    #cl = ax1.getp(c, 'xmajorticklabels')
    c.ax.tick_params(labelsize=20)

    pad = 15  # in points

    ax1.annotate(pltlbl, xy=(0, 0.5), xytext=(-ax1.yaxis.labelpad - pad, 0),
                 fontsize=32, fontname='Arial',
                 xycoords=ax1.yaxis.label, textcoords='offset points',
                 ha='right', va='center', rotation=90)

    plt.show()
    plt.subplots_adjust(hspace=0, wspace=0.05, top=.97, bottom=0.15, left=.05,
                        right=.95)


def latweightmean(data, msk='no'):
    from netcdfread import ncread
    lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')[0, 0, :]
    lsm1 = np.ones((145, 192))
    for i in range(145):
        for j in range(192):
            if lsm[i, j] == 0:
                lsm1[i, j] = 0
    lon = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'longitude0')
    lat = ncread('/network/aopp/hera/mad/bakerh/HAPPI/batch_518/atmos/item3236_monthly_mean/item3236_monthly_mean_a011_2006-01_2016-12.nc', 'latitude0')
    meshlat = np.zeros([np.ma.size(lon), np.ma.size(lat)])
    meshlat[:, :] = lat
    meshlatweight = np.transpose(np.cos(meshlat * np.pi/180))
    weighted = data * meshlatweight
    if msk == 'yes':
        mask = np.zeros((145, 192))
        mask[:] = lsm1
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'NH':
        mask = np.zeros((73, 192))
        lsm1 = lsm1[:73, :]
        mask[:] = lsm1
        weighted = weighted[:73, :]
        meshlatweight = meshlatweight[:73, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'TROP':
        mask = np.zeros((49, 192))
        lsm1 = lsm1[48:97, :]
        mask[:] = lsm1
        weighted = weighted[48:97, :]
        meshlatweight = meshlatweight[48:97, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'MID':
        mask = np.zeros((48, 192))
        lsm1 = lsm1[:48, :]
        mask[:] = lsm1
        weighted = weighted[:48, :]
        meshlatweight = meshlatweight[:48, :]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    elif msk == 'no':
        meaned = np.mean(weighted) / np.mean(meshlatweight)
    elif np.ma.size(msk) == 5:
        mask = np.zeros((msk[2]-msk[1],
                         192+msk[4]-msk[3]))
        lsm1 = np.concatenate((lsm1[msk[1]:msk[2], msk[3]:], lsm1[msk[1]:msk[2], :msk[4]]), axis=1)
        mask[:] = lsm1
        weighted = np.concatenate((weighted[msk[1]:msk[2], msk[3]:], weighted[:, msk[1]:msk[2], :msk[4]]), axis=2)
        meshlatweight = np.concatenate((meshlatweight[msk[1]:msk[2], msk[3]:], meshlatweight[msk[1]:msk[2], :msk[4]]), axis=1)
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    else:
        mask = np.zeros((msk[1]-msk[0], msk[3]-msk[2]))
        lsm1 = lsm1[msk[0]:msk[1], msk[2]:msk[3]]
        mask[:] = lsm1
        weighted = weighted[msk[0]:msk[1], msk[2]:msk[3]]
        meshlatweight = meshlatweight[msk[0]:msk[1], msk[2]:msk[3]]
        meaned = np.mean(weighted[mask > 0.5]) / np.mean(meshlatweight[lsm1 > 0.5])
    return meaned


def spreadplt(tx90p, tcres, global_tas):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.colors import Normalize
    import matplotlib.cm as cm
    plt.figure()
    time = np.arange(2006, 2101, 1)
    cmap = cm.viridis
    norm = Normalize(vmin=1.1, vmax=np.max(tcres))
    gt_mov = np.zeros((13, 86))
    tx90p_mov = np.zeros((13, 86))
    for i in range(86):
        gt_mov[:, i] = np.mean(global_tas[:, 144+i:i+154], axis=1)
        tx90p_mov[:, i] = np.mean(tx90p[:, i:i+10], axis=1)
    for i in range(13):
        if tcres[i] != 0:
            plt.plot(gt_mov[i], tx90p_mov[i], color=cmap(norm(tcres[i])))
    #plt.xlim([2006, 2100])
    plt.ylim([10, 90])
    #plt.xlabel('Year')
    plt.ylabel('TX90p (days)')




















