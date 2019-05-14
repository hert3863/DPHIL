# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 14:34:10 2016

@author: bakerh

NAME
    GCM Script Writer
PURPOSE
    Creates scripts for GCM run with varying parameters
"""


def alldaily(experiment, run, invariable):
    '''
    Takes an experiment name as an input
    and returns the full daily data
    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variable : array
        Array of variable specified
    '''
    def ncreadfield(filelocation, invariable):
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

    import numpy as np
    import os
    npu = np.zeros([1, 37, 64])
    for d in os.listdir('/network/aopp/hera/mad/bakerh/fms_tmp/' +
                        experiment + '/' + run + '/combine/'):
        u = ncreadfield('/network/aopp/hera/mad/bakerh/fms_tmp/' +
                        experiment + '/' + run + '/combine/' + d, invariable)
        u = np.mean(u, axis=3)
        npu = np.concatenate((npu, u), axis=0)
        del u
    npu = npu[1:, :, :]
    variable = np.zeros([int(np.ma.size(npu, axis=0)/4), 37, 64])
    for i in range(np.ma.size(variable, axis=0)):
        variable[i, :] = (npu[4*i, :] + npu[4*i+1, :] +
                          npu[4*i+2, :] + npu[4*i+3, :])/4
    del npu
    return variable


def combine1(experiment, run, invariable):
    '''
    outputs numpy array combining netcdfs for a run

    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    import numpy as np
    import os
    from netCDF4 import Dataset
    variable = np.zeros((37, 64))
    s = len(os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                       experiment + '/' + run + '/history/')) - 1
    for d in os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        experiment + '/' + run + '/history/'):
        if d.find('daily.nc') == -1:
            nc_f = ('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                    experiment + '/' + run + '/history/' + d)
            nc_fid = Dataset(nc_f, 'r')
            data = nc_fid.variables[invariable][:]
            variable = data + variable
    variable = variable / s
    return variable


def dailyanalysis(experiment):
    """
    Takes an experiment name as an input
    and completes the daily analysis

    Parameters
    ----------
    experimentname: str
        name of experiment
    """
    import os
    for fn in os.listdir('/network/aopp/hera/mad/bakerh/fms_tmp/' +
                         experiment):
        if fn.find('exe.fms') == -1 and fn.find('mppnccombine.ifc') == -1:
            storedaily('/network/aopp/hera/mad/bakerh/fms_tmp/' + experiment +
                       '/' + fn + '/combine/',
                       '/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                       experiment + '/' + fn + '/history/')
            print('Completed ' + fn)


def generate(experiment, ifilename, parameterarray):
    """
    Takes a default input file, and an array of parameters.
    The output file names are named the same as the runname
    to avoid confusion. A final script is then generated
    which contains all the output filenames ready to execute

    Parameters
    ----------
    experiment: str
        name of experiment
    ifilename: str
        name of input file
    parameterarray: array
        parameters for run with runnames on each row, parameters
        on columns
    """
    import numpy as np
    import os
    # create file in fms_tmp and copy in requisite files
    rsyncstring = "rsync -a --exclude='climspinup' \
'/network/aopp/hera/mad/bakerh/fms_tmp/climspinup/' \
'/network/aopp/hera/mad/bakerh/fms_tmp/" + experiment + "'"
    os.system(rsyncstring)
    # separate code to change run_names and write initial files
    runfile = open('/home/bakerh/fms/exp/' + experiment +
                   '/run/' + 'runfile', 'w')
    runfile.write('#!/bin/csh -f\n')
    for i in range(np.ma.size(parameterarray, axis=0)-1):
        ifile = open('/home/bakerh/fms/exp/' + experiment +
                     '/run/' + ifilename, 'r')
        lines = ifile.readlines()
        ifile.close()
        ofile = open('/home/bakerh/fms/exp/' + experiment + '/run/' +
                     parameterarray[i+1, 0], 'w')
        for line in lines:
            if line.find('label for') != -1:
                ofile.write('set run_name = ' + parameterarray[i+1, 0] + '\n')
            else:
                ofile.write(line)
        ofile.close()
        os.chmod('/home/bakerh/fms/exp/' + experiment + '/run/' +
                 parameterarray[i+1, 0], 33279)
        runfile.write('./' + parameterarray[i+1, 0] + '\n')
        # copy restart file and create restart text file
        dirtomake = "mkdir '/network/aopp/hera/mad/bakerh/fms_tmp/\
" + experiment + "/" + parameterarray[i+1, 0] + "'"
        os.system(dirtomake)
        copyrestart = "rsync -a '/network/aopp/hera/mad/bakerh/fms_tmp/\
climspinup/climspinup/output/restart/day3600h00.cpio' \
'/network/aopp/hera/mad/bakerh/fms_tmp/\
" + experiment + "/" + parameterarray[i+1, 0] + "'"
        os.system(copyrestart)
        rfile = open('/network/aopp/hera/mad/bakerh/fms_tmp/' + experiment +
                     '/' + parameterarray[i+1, 0] + '/reload_commands', 'w')
        rfile.write('set irun         =  1\n\
set init_cond    =  /network/aopp/hera/mad/bakerh/fms_tmp/' +
                    experiment + '/' + parameterarray[i+1, 0] +
                    '/day3600h00.cpio \nset ireload      =  2')
        rfile.close()
    runfile.close()
    os.chmod('/home/bakerh/fms/exp/' + experiment + '/run/' + 'runfile', 33279)
    # now alter parameters
    for i in range(np.ma.size(parameterarray, axis=0)-1):
        for j in range(np.ma.size(parameterarray, axis=1)-1):
            parameters('/home/bakerh/fms/exp/' + experiment +
                       '/run/' + parameterarray[i+1, 0],
                       '/home/bakerh/fms/exp/' +
                       experiment + '/run/' + parameterarray[i+1, 0],
                       parameterarray[0, j+1], parameterarray[i+1, j+1])


def indiv1(experiment, run, invariable):
    '''
    outputs numpy array combining netcdfs for a run in slices

    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    import numpy as np
    import os
    from netCDF4 import Dataset
    yrs = []
    for d in os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                        experiment + '/run' + str(run) + '/history/'):
        if d.find('daily.nc') == -1 and d.find('.', 0, 1) == -1:
            yrs.append(d)
    variable = np.zeros((37, 64))
    variable = np.zeros((len(yrs), 37, 64))
    for i in range(len(yrs)):
        nc_f = ('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                experiment + '/run' + str(run) + '/history/' + yrs[i])
        nc_fid = Dataset(nc_f, 'r')
        data = nc_fid.variables[invariable][:]
        variable[i, :, :] = data
    return variable


def multiread1(experiment, invariable):
    '''
    outputs numpy array of all runs in experiment

    Parameters
    ----------
    experiment : str
        experimentname
    run: int
        run number
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    import numpy as np
    import os
    s = len(os.listdir('/network/aopp/hera/mad/bakerh/data/FMS/output/' +
                       experiment))
    variables = np.zeros((s, 37, 64))
    for i in range(len(os.listdir('/network/aopp/hera/mad/' +
                                  'bakerh/data/FMS/output/' + experiment))):
        variables[i, :, :] = combine1(experiment, 'run' + str(i), invariable)
    return variables


def parameterarray():
    """
    Produces the parameter array by taking using
    specified inputs

    Returns
    ----------
    parameterarray: array
        parameters for run with runnames on each row, parameters
        on columns
    """
    import numpy as np
    numberofruns = input('Enter number of runs: ')
    numberofpara = input('Enter number of varying parameters: ')
    paraarray = np.full([int(numberofruns)+1, int(numberofpara)+1], 0)
    paraarray = np.ndarray.astype(paraarray, str)
    paraarray[0, 0] = 'runname'
    for i in range(int(numberofruns)):
        paraarray[i+1, 0] = 'run' + str(i+1)
    for j in range(int(numberofpara)):
        parameter = input('Enter parameter ' + str(j+1) + ': ')
        paraarray[0, j+1] = parameter
    # if paraarray[0, 1] == 'x_0' and paraarray[0, 2] == 'y_0':
    return paraarray


def parameters(ifilename, ofilename, parametername, parameter):
    """
    Alters parameters in script

    Parameters
    ----------
    ifilename: str
        name of input file
    ofilename: str
        output file name
    parametername: str
        parameter name to alter
    parameter: array
        parameter for run
    """
    ifile = open(ifilename, 'r')
    lines = ifile.readlines()
    ifile.close()
    ofile = open(ofilename, 'w')
    for line in lines:
        if line.find(parametername) != -1:
            ofile.write('\t' + parametername + '=' + str(parameter) + ',\n')
        else:
            ofile.write(line)
    ofile.close()


def storedaily(tempdir, outputdir):
    '''
    takes the temporary datafiles from the GCM runs
    and takes out the u, v and T data, turns it into daily
    data and exports to netCDF at 3 sigma levels

    Parameters
    ----------
    tempdir : str
        directory of temporary datafiles
    '''
    import os
    import numpy as np

    def ncreadfield(filelocation, invariable):
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
        variable = nc_fid.variables[invariable][:, [18, 24, 31], :, :]
        return variable

    from netCDF4 import Dataset

    ofile = Dataset(outputdir + 'daily.nc',
                    'w', format='NETCDF4')
    level = ofile.createDimension('level', 3)
    time = ofile.createDimension('time', None)
    lat = ofile.createDimension('lat', 64)
    times = ofile.createVariable('time', 'f8', ('time',))
    levels = ofile.createVariable('level', 'f4', ('level',))
    latitudes = ofile.createVariable('lat', 'f4', ('lat',))
    uwnd = ofile.createVariable('uwnd', 'f4', ('time', 'level', 'lat',))
    vwnd = ofile.createVariable('vwnd', 'f4', ('time', 'level', 'lat',))
    temp = ofile.createVariable('temp', 'f4', ('time', 'level', 'lat',))
    nptimes = np.zeros(0)
    for fn in os.listdir(tempdir):
        nc_f = tempdir + str(fn)  # filename
        nc_fid = Dataset(nc_f, 'r')
        nplevels = nc_fid.variables['pfull'][[18, 24, 31]]  # for 37 sigma
        nplat = nc_fid.variables['lat'][:]
        nptime = nc_fid.variables['time'][:]
        nptimes = np.concatenate((nptimes, nptime), axis=0)
    nptimes = nptimes[1:]
    nptimesdaily = nptimes[::4]
    times[:] = nptimesdaily
    latitudes[:] = nplat
    levels[:] = nplevels

    npu = np.zeros([1, 3, 64])
    for fn in os.listdir(tempdir):
        u = ncreadfield(tempdir + str(fn), 'ucomp')
        u = np.mean(u, axis=3)
        npu = np.concatenate((npu, u), axis=0)
        del u
    npu = npu[1:, :, :]
    npudaily = np.zeros([int(np.ma.size(npu, axis=0)/4), 3, 64])
    for i in range(np.ma.size(npudaily, axis=0)):
        npudaily[i, :] = (npu[4*i, :] + npu[4*i+1, :] +
                          npu[4*i+2, :] + npu[4*i+3, :])/4
    uwnd[:] = npudaily
    del npu
    del npudaily

    npv = np.zeros([1, 3, 64])
    for fn in os.listdir(tempdir):
        v = ncreadfield(tempdir + str(fn), 'vcomp')
        v = np.mean(v, axis=3)
        npv = np.concatenate((npv, v), axis=0)
        del v
    npv = npv[1:, :, :]
    npvdaily = np.zeros([int(np.ma.size(npv, axis=0)/4), 3, 64])
    for i in range(np.ma.size(npvdaily, axis=0)):
        npvdaily[i, :] = (npv[4*i, :] + npv[4*i+1, :] +
                          npv[4*i+2, :] + npv[4*i+3, :])/4
    vwnd[:] = npvdaily
    del npv
    del npvdaily

    npt = np.zeros([1, 3, 64])
    for fn in os.listdir(tempdir):
        t = ncreadfield(tempdir + str(fn), 'temp')
        t = np.mean(t, axis=3)
        npt = np.concatenate((npt, t), axis=0)
        del t
    npt = npt[1:, :, :]
    nptdaily = np.zeros([int(np.ma.size(npt, axis=0)/4), 3, 64])
    for i in range(np.ma.size(nptdaily, axis=0)):
        nptdaily[i, :] = (npt[4*i, :] + npt[4*i+1, :] +
                          npt[4*i+2, :] + npt[4*i+3, :])/4
    temp[:] = nptdaily
    del npt
    del nptdaily

    ofile.close()
