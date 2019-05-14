# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:34:02 2016

@author: bakerh
"""
import numpy as np


def mapplot(plotdata, lat, lon, title=''):
    """
    Plots input grid with map of world overlaid

    Parameters
    ----------
    plotdata: array
        data being plotted
    title: str
        optional title
    """
    from mpl_toolkits.basemap import Basemap, shiftgrid
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    '''
    long = np.zeros((193))
    long[:-1], long[-1] = lon, 360.000
    temp = plotdata[:, 0]
    plotdata = np.c_[plotdata, temp]
    '''
    plotdata, lon = shiftgrid(180., plotdata, lon, start=False)
    meshlon, meshlat = np.meshgrid(lon, lat)
    plt.figure()
    m = Basemap(projection='cyl', llcrnrlat=-90, urcrnrlat=90,
                llcrnrlon=-180, urcrnrlon=180, resolution='c')
    m.drawcoastlines()
    m.drawmapboundary()
    x, y = m(meshlon, meshlat)
    my_cmap = plt.cm.jet(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    # ctrs = np.linspace(-30, 30, 21)
    ctrs = np.linspace(-20, 50, 29)
    # ctrs = np.linspace(-0.2, 1, 12)
    for i in range(np.ma.size(lat)):
        for j in range(np.ma.size(lon)):
            if plotdata[i, j] == 'masked':
                plotdata[i, j] = np.nan
    plot = m.contourf(x, y, plotdata, ctrs,
                      cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                      extend='both')
    plt.colorbar(plot, orientation='horizontal',
                 shrink=0.5, spacing='proportional')
    parallels = m.drawparallels(np.arange(-90., 91., 15.))
    meridians = m.drawmeridians(np.arange(-180., 181., 30.))
    m.drawparallels(parallels, labels=[True, True, True, True])
    m.drawmeridians(meridians, labels=[True, True, True, True])
    plt.title(title, y=1.08)
    plt.show()


def colourscale(plotdata):
    """
    Takes data being plotted and normalises the colourscale between largest
    data value and its negative multiple

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


def ancil(inputsst, inputsice, outputfile):
    '''
    takes the inputnpy and exports to ancillary

    Parameters
    ----------
    inputsst: sst array
        sst data to be exported
    inputsice: sice array
        sice data to be exported
    outputfile: str
        directory and filename for output
    '''

    def npy2netcdf(inputsst, inputsice, outputfile):
        '''
        takes the inputnpy and exports to netCDF

        Parameters
        ----------
        inputsst: sst array
            sst data to be exported
        inputsice: sice array
            sice data to be exported
        outputfile: str
            directory and filename for output
        '''
        import numpy as np
        from netCDF4 import Dataset

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
            from netCDF4 import Dataset
            nc_f = filelocation  # filename
            nc_fid = Dataset(nc_f, 'r')
            variable = nc_fid.variables[invariable][:]
            return variable
        # invert latitudes
        inputsst = inputsst[:, ::-1, :]
        inputsice = inputsice[:, ::-1, :]
        # convert to 14 months
        inputsst = np.concatenate((np.expand_dims(inputsst[11, :], 0),
                                   inputsst,
                                   np.expand_dims(inputsst[1, :], 0)))
        inputsice = np.concatenate((np.expand_dims(inputsice[11, :], 0),
                                    inputsice,
                                    np.expand_dims(inputsice[1, :], 0)))
        inputsst = inputsst + 273.15
        dates = np.linspace(-15, 765, 27)
        dates5 = np.linspace(-12.5, 762.5, 156)
        lons = np.linspace(0, 358.125, 192)
        lats = np.linspace(90, -90, 145)
        # interpolate
        outputsst = np.zeros((np.ma.size(inputsst, 0)*6-6, 145, 192))
        outputsice = np.zeros((np.ma.size(inputsice, 0)*6-6, 145, 192))
        for i in range(145):
            for j in range(192):
                outputsst[:, i, j] = np.interp(dates5, dates,
                                               inputsst[:, i, j])
                outputsice[:, i, j] = np.interp(dates5, dates,
                                                inputsice[:, i, j])
        outputsst = outputsst[3:, :]
        outputsice = outputsice[3:, :]
        # apply landsea mask
        lsm = ncread('/home/bakerh/Documents/DPhil/CPDN/\
Weather-at-Home_ancilmaker-master/lsm_n96_add.nc', 'lsm')
        lsm = lsm[0, 0, :]
        lsm1 = np.zeros((np.ma.size(outputsst, axis=0), 145, 192))
        for i in range(np.ma.size(outputsst, axis=0)):
            lsm1[i, :] = lsm
        masksst = np.ma.masked_array(outputsst, lsm1)    
        masksice = np.ma.masked_array(outputsice, lsm1)
        # begin saving
        ofile = Dataset('/network/aopp/hera/mad/bakerh/RPM/SSTancils/' +
                        outputfile + '.nc', 'w', format='NETCDF3_CLASSIC')
        ofile.createDimension('time', None)
        ofile.createDimension('lat', 145)
        ofile.createDimension('lon', 192)
        times = ofile.createVariable('time', 'f4', ('time',))
        times.units = 'days since 1998-12-01 00:00:00'
        times.calendar = '360_day'
        longitudes = ofile.createVariable('lon', 'f4', ('lon',))
        longitudes.units = 'degrees_east'
        latitudes = ofile.createVariable('lat', 'f4', ('lat',))
        latitudes.units = 'degrees_north'
        ofile.createDimension('surface', 1)
        surface = ofile.createVariable('surface', 'f4', ('surface',))
        surface.units = 'level'
        sst = ofile.createVariable('sst', 'f4', ('time', 'surface', 'lat', 'lon'),
                                   fill_value=2e20)
        sst.units = 'K'
        sice = ofile.createVariable('sice', 'f4', ('time', 'surface', 'lat', 'lon'),
                            fill_value=2e20)
        sice.units = 'fraction'
        surface[:] = 0
        latitudes[:] = lats
        longitudes[:] = lons
        times[:] = dates5[3:]
        sst[:] = masksst
        sice[:] = masksice

        ofile.close()

    def netcdf2ancil():
        '''
        converts netcdf3 to ancil
        using Mitch's code
        '''
        import os
        cline = ('python /home/bakerh/Documents\
/DPhil/Python/CPDN/gen_cpdn_sst_sice_ancil.py /network/aopp/hera/mad/bakerh/RPM/SSTancils/' + outputfile + '.nc sst /network/aopp/hera/mad/bakerh/RPM/SSTancils/' + outputfile + '.nc sice 1998')
        os.system(cline)
        cpy = "mv '/home/bakerh/Documents\
/DPhil/Python/sst_ancil' '/network/aopp/hera/mad/bakerh/RPM/SSTancils/" + outputfile + "'"
        os.system(cpy)
        cpy1 = "mv '/home/bakerh/Documents\
/DPhil/Python/sice_ancil' '/network/aopp/hera/mad/bakerh/RPM/SSTancils/sice_ancil'"
        os.system(cpy1)
        rem = "rm '/network/aopp/hera/mad/bakerh/RPM/SSTancils/" + outputfile + ".nc'"
        os.system(rem)
        gzip = "gzip '/network/aopp/hera/mad/bakerh/RPM/SSTancils/" + outputfile + "'"
        os.system(gzip)

    npy2netcdf(inputsst, inputsice, outputfile)
    netcdf2ancil()

'''
for i in range(500,10001):
    ancil(ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert'+str(i)+'_c040926.nc','SST'),ncread('/network/aopp/hera/mad/bakerh/RPM/SSTpatches/sst_HadOIBl_bc_N96_clim_pert'+str(i)+'_c040926.nc','ice_cov'),'sst_HadOIBl_bc_N96_clim_pert'+str(i))
'''


def genpertlist(include_zero_pert=False):
    # generate a list of possible perturbation files
    import random
    pert_list = []
    pert_base_year = "ic1961"
    pert_suffix = "_N96"
    scale_set = ["10", "11", "12", "14", "16"]

    for mm in range(1, 13):
        for dd in range(1, 30):
            for sc in range(0, 5):
                pert_str = pert_base_year + "%02d" % mm + "%02d_" % dd + \
                              scale_set[sc] + pert_suffix
                pert_list.append(pert_str)

    # shuffle the list so it has a random order of perturbations
    random.shuffle(pert_list)
    # add the zero perturbation to the front of the list
    if include_zero_pert:
        pert_list.insert(0, "ic00000000_10_N96")

    return pert_list


def genxmlwu(ensemble_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('b001')
    pert_list = genpertlist()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwu', 'w')
    runfilecont = open('/home/bakerh/Documents/DPhil/CPDN/xmlwucontrol', 'w')
    for i in range(1, ensemble_size+1):
        runfilecont.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>restart_atmos_rpm_5yr</file_atmos>\n\
                <file_sst>sst_HadOIBl_bc_N96_clim</file_sst>\n\
                <file_sice>sice_HadOIBl_bc_N96_clim</file_sice>\n\
                <file_so2dms>so2dms_1980_2000_av</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_1980_2000_av</file_ozone>\n\
                <file_pert>' + pert_list[np.remainder(i,1740).astype(int)] + '</file_pert>\n\
                <file_solar>solar_1980_2000_av</file_solar>\n\
                <file_volcanic>volc_allzero</file_volcanic>\n\
                <file_ghg>ghg_1980_2000_av</file_ghg>\n\
                <run_years>2</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_year>1998</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>restart_atmos_rpm_5yr</file_atmos>\n\
                <file_sst>sst_HadOIBl_bc_N96_clim_pert' + str(i) + '</file_sst>\n\
                <file_sice>sice_HadOIBl_bc_N96_clim</file_sice>\n\
                <file_so2dms>so2dms_1980_2000_av</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_1980_2000_av</file_ozone>\n\
                <file_pert>' + pert_list[np.remainder(i,1740).astype(int)] + '</file_pert>\n\
                <file_solar>solar_1980_2000_av</file_solar>\n\
                <file_volcanic>volc_allzero</file_volcanic>\n\
                <file_ghg>ghg_1980_2000_av</file_ghg>\n\
                <run_years>2</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_year>1998</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('RP' + str(i), pert_list[np.remainder(i,1740).astype(int)])
        anc.Next()
    runfile.close()
    runfilecont.close()
    return wuinfo
