# -*- coding: utf-8 -*-
'''
Created on Wed Feb 10 11:25:04 2016

@author: bakerh

NAME
    NetCDF with Python
PURPOSE
    To read data with NetCDF files using
    a NetCDF file from the Reanalyses.
'''


def ncdump(filelocation, verb=True):
    '''
    ncdump outputs dimensions, variables and their attribute information.
    The information is similar to that of NCAR's ncdump utility.
    ncdump requires a valid instance of Dataset.

    Parameters
    ----------
    fileloaction : str
    netCDF4.Dataset file location
    verb : Boolean
        whether or not nc_attrs, nc_dims, nc_vars and nc_lst are printed

    Returns
    -------
    nc_attrs : list
        A Python list of the NetCDF file global attributes
    nc_dims : list
        A Python list of the NetCDF file dimensions
    nc_vars : list
        A Python list of the NetCDF file variables
    nc_lst : list
        List of all variables with out details
    '''
    def print_ncattr(key):
        """
        Prints the NetCDF file attributes for a given key

        Parameters
        ----------
        key : unicode
            a valid netCDF4.Dataset.variables key
        """
        try:
            print("\t\ttype:", repr(nc_fid.variables[key].dtype))
            for ncattr in nc_fid.variables[key].ncattrs():
                print('\t\t%s:' % ncattr,
                      repr(nc_fid.variables[key].getncattr(ncattr)))
        except KeyError:
            print("\t\tWARNING: %s does not contain variable attributes" % key)
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    nc_f = filelocation  # filename
    nc_fid = Dataset(nc_f, 'r')
    # Dataset is the class behavior to open
    # and create an instance of the ncCDF4 class
    # NetCDF global attributes
    nc_attrs = nc_fid.ncattrs()
    if verb:
        print("NetCDF Global Attributes:")
        for nc_attr in nc_attrs:
            print('\t%s:' % nc_attr, repr(nc_fid.getncattr(nc_attr)))
    nc_dims = [dim for dim in nc_fid.dimensions]  # list of nc dimensions
    # Dimension shape information.
    if verb:
        print("NetCDF dimension information:")
        for dim in nc_dims:
            print("\tName:", dim)
            print("\t\tsize:", len(nc_fid.dimensions[dim]))
            print_ncattr(dim)
    # Variable information.
    nc_vars = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable information:")
        for var in nc_vars:
            if var not in nc_dims:
                print('\tName:', var)
                print("\t\tdimensions:", nc_fid.variables[var].dimensions)
                print("\t\tsize:", nc_fid.variables[var].size)
                print_ncattr(var)
    # Variable list
    '''
    nc_lst = [var for var in nc_fid.variables]  # list of nc variables
    if verb:
        print("NetCDF variable list:")
        for var in nc_lst:
            print('\tName:', var, '  Long Name:',
                  nc_fid.variables[var].long_name)
    '''


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


def ncsave(path, t, latitude, longitude, data, fldname):
    #ncsave(path, latitude, longitude, hgts, t, data, reftime, fldname):
    '''
    reftime is sting in format yyyy-mm-dd
    '''
    from netCDF4 import Dataset
    ofile = Dataset(path + '.nc', 'w', format='NETCDF3_64BIT_OFFSET')
    ofile.createDimension('time', None)
    ofile.createDimension('lat', len(latitude))
    ofile.createDimension('lon', len(longitude))
    #ofile.createDimension('level', len(hgts))
    times = ofile.createVariable('time', 'f8', ('time',))
    #times.units = 'days since ' + reftime + ' 00:00:00'
    #times.calendar = '360_day'
    longitudes = ofile.createVariable('lon', 'f4', ('lon',))
    longitudes.units = 'degrees_east'
    latitudes = ofile.createVariable('lat', 'f4', ('lat',))
    latitudes.units = 'degrees_north'
    #levels = ofile.createVariable('level', 'f4', ('level',))
    #levels.units = 'hPa'
    field = ofile.createVariable(fldname, 'f4', ('time', 'lat',
                                                 'lon'),
                                 fill_value=2e20)
    latitudes[:] = latitude
    longitudes[:] = longitude
    times[:] = t
    #levels[:] = hgts
    field[:] = data

    ofile.close()


def ncreadmulti(filelocation, invariable, start, end):
    '''
    ncread outputs numpy arrays of called variables in netCDF4 file

    Parameters
    ----------
    fileloaction : str
        NetCDFfile directory
    invariable : str
        variable stored in file to be imported
    start: int
        start file to import
    end: int
        end file to import

    Returns
    -------
    variable : array
        Array of variable specified
    '''
    import numpy as np
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    variable = np.zeros((end-start+1, 37, 64))
    for i in range(end-start+1):
        nc_f = filelocation + '/day0' + repr(i+start) + 'h00.nc'
        nc_fid = Dataset(nc_f, 'r')
        variable[i, :, :] = nc_fid.variables[invariable][:]
    return variable
