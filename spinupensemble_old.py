# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 10:02:07 2016

@author: bakerh

NAME
    Spin up ensemble
PURPOSE
    All functions necessary to setup and analyse
    the spin up ensemble experiment.
"""
import numpy as np


def restartstates():
    '''
    Takes an initial default file with sigma and lat coordinate and
    creates the runfiles and restart files from
    the spinupstates experiment
    '''
    import os
    dirtomak = "mkdir '/network/aopp/hera/mad/bakerh/fms_tmp/spinupensembleej2/'"
    os.system(dirtomak)
    # code to change run_names and write initial files
    runfile = open('/home/bakerh/fms/exp/spinupensembleej2/run/runfile', 'w')
    runfile.write('#!/bin/csh -f\n')
    for i in range(200):
        ifile = open('/home/bakerh/fms/exp/spinupensembleej2/run/default', 'r')
        lines = ifile.readlines()
        ifile.close()
        ofile = open('/home/bakerh/fms/exp/spinupensembleej2/run/member' +
                     str(i+1), 'w')
        for line in lines:
            if line.find('label for') != -1:
                ofile.write('set run_name = member' + str(i+1) + '\n')
            else:
                ofile.write(line)
        ofile.close()
        os.chmod('/home/bakerh/fms/exp/spinupensembleej2/run/member' +
                 str(i+1), 33279)
        runfile.write('./member' + str(i+1) + '\n')
        # copy restart file and create restart text file
        dirtomake = "mkdir '/network/aopp/hera/mad/bakerh/fms_tmp/\
spinupensembleej2/member" + str(i+1) + "'"
        os.system(dirtomake)
        copyrestrt = "rsync -a --no-R --no-implied-dirs '/network/aopp/hera/mad/bakerh/fms_tmp/\
spinupstates/spinupstates/output/restart/day" + str(3610+10*i) + "h00.cpiotmp' \
'/network/aopp/hera/mad/bakerh/fms_tmp/\
spinupensembleej2/member" + str(i+1) + "/day" + str(3610+10*i) + "h00.cpio'"
        os.system(copyrestrt)
        rfile = open('/network/aopp/hera/mad/bakerh/fms_tmp/spinupensembleej2/' +
                     'member' + str(i+1) + '/reload_commands', 'w')
        rfile.write('set irun         =  1\n\
set init_cond    =  /network/aopp/hera/mad/bakerh/fms_tmp/spinupensembleej2/\
member' + str(i+1) + '/day' + str(3610+10*i) + 'h00.cpio \nset ireload   =  2')
        rfile.close()
    runfile.close()
    os.chmod('/home/bakerh/fms/exp/spinupensembleej2/run/runfile', 33279)


def controlimport(invariable):
    '''
    imports and means over the control run

    Parameters
    ----------
    invariable : str
        variable stored in file to be imported

    Returns
    -------
    variables : array
        Array of variable specified
    '''
    import os
    from netCDF4 import Dataset
    variable = np.zeros((64))
    for d in os.listdir('/network/aopp/hera/mad/bakerh/fms_tmp/spinupstates/' +
                        'spinupstates/combine/'):
        nc_f = ('/network/aopp/hera/mad/bakerh/fms_tmp/spinupstates/' +
                'spinupstates/combine/' + d)
        nc_fid = Dataset(nc_f, 'r')
        data = nc_fid.variables[invariable][:]
        data = np.mean(data, axis=2)
        data = np.mean(data, axis=0)
        variable = data + variable
    variable = variable / 200
    return variable


def ensembleimport(exp):
    '''
    Imports ensemble quantities
    '''
    import glob
    from netCDF4 import Dataset
    u = np.zeros((720, 37, 64))
    v = np.zeros((720, 37, 64))
    temp = np.zeros((720, 37, 64))
    d = glob.glob('/network/aopp/hera/mad/bakerh/fms_tmp/' + exp +
                  '/member*')
    for e, f in enumerate(d):
        nc_f = glob.glob(f + '/combine/*')
        nc_fid = Dataset(nc_f[0], 'r')
        data1 = np.mean(nc_fid.variables['ucomp'][:], axis=3)
        data2 = np.mean(nc_fid.variables['vcomp'][:], axis=3)
        data3 = np.mean(nc_fid.variables['temp'][:], axis=3)
        u = data1 + u
        v = data2 + v
        temp = data3 + temp
        print(str(e))
    u = u / 200
    v = v / 200
    temp = temp / 200
    return u, v, temp


def ptemp(temp, pfull):
    '''
    computes the potential temperature
    '''
    k = 0.286
    pottemp = np.zeros((37, 64))
    for a in range(64):
        pottemp[:, a] = temp[:, a] * (1000 / pfull) ** k
    return pottemp


def sfctn(v, phalf, lat):
    '''
    computes the streamfunction
    '''
    c = -200 * np.pi * 6371000 / 9.81
    streamfctn = np.zeros((37, 64))
    meshlat, meshp = np.meshgrid(lat, phalf[:-1])
    vw = v * np.cos(np.pi*meshlat/180)
    for i in range(37):
        streamfctn[i, :] = streamfctn[i-1, :] + c * (vw[i, :] *
                                                     (phalf[i+1]-phalf[i]))
    return streamfctn


def geoz(temp, phalf, pfull):
    '''
    computes geopotential height
    '''
    zhalf = np.zeros((38, 64))
    z = np.zeros((37, 64))
    ln_phalf = np.log(phalf*100)
    ln_phalf[0] = 0.0
    ln_pfull = np.log(pfull*100)
    for i in range(36, 0, -1):
        zhalf[i, :] = zhalf[i+1, :] + 287.04 * temp[i, :] * (ln_phalf[i+1] -
                                                             ln_phalf[i]) / 9.8
    for i in range(37):
        z[i, :] = zhalf[i+1, :] + 287.04 * temp[i, :] * (ln_phalf[i+1] -
                                                         ln_pfull[i]) / 9.8
    return z


def lapserate(t, z, sigma, lat):
    """
    Produces plot of lapse rate of T data
    """
    import numpy as np
    dT = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    dz = np.zeros((np.ma.size(sigma), np.ma.size(lat)))
    for i in range(np.ma.size(sigma, axis=0)-1):
        dT[i, :] = t[i+1, :] - t[i, :]
    for i in range(np.ma.size(sigma, axis=0)-1):
        dz[i, :] = z[i+1, :] - z[i, :]
    lapse = -1000 * dT[0:-1] / dz[0:-1]
    return lapse


def prime(exp):
    '''
    Computes the transients
    NEED COSINE WEIGHTING
    '''
    import glob
    from netCDF4 import Dataset  # http://code.google.com/p/netcdf4-python/
    uv = np.zeros((720, 37, 64))
    vT = np.zeros((720, 37, 64))
    for i in range(720):
        u = np.zeros((200, 37, 64, 128))
        v = np.zeros((200, 37, 64, 128))
        T = np.zeros((200, 37, 64, 128))
        d = glob.glob('/network/aopp/hera/mad/bakerh/fms_tmp/' + exp +
                      '/member*')
        for e, f in enumerate(d):
            nc_f = glob.glob(f + '/combine/*')
            nc_fid = Dataset(nc_f[0], 'r')
            u[e] = nc_fid.variables['ucomp'][i, :]
            v[e] = nc_fid.variables['vcomp'][i, :]
            T[e] = nc_fid.variables['temp'][i, :]
        uv1 = (u - np.mean(u, axis=0)) * (v - np.mean(v, axis=0))
        vT1 = (T - np.mean(T, axis=0)) * (v - np.mean(v, axis=0))
        uv1 = np.mean(uv1, axis=3)
        vT1 = np.mean(vT1, axis=3)
        uv1 = np.mean(uv1, axis=0)
        vT1 = np.mean(vT1, axis=0)
        uv[i, :] = uv1
        vT[i, :] = vT1
        print(str(i))
    return uv, vT


def diverg(field, lat, pfull):
    meshlat, meshp = np.meshgrid(lat, pfull)
    phi = meshlat*np.pi/180
    cs_phi = np.cos(phi)
    dphi = lat*np.pi/180
    a = 6.371e6
    div = np.zeros(np.shape(field))
    for i in range(len(field)):
        div[i] = np.gradient(field[i]*cs_phi, dphi, axis=1) / (a*cs_phi)
    return div


def lanczos(field, wi=51, co=0.0125, dim=1):
    def lweights(window, cutoff):
        """
        Calculate weights for a low pass Lanczos filter
        and then applies to time series

        Parameters
        ----------

        window: int
            The length of the filter window.
        cutoff: float
            The cutoff frequency in inverse time steps.

        Returns
        -------
        w: array
            filter weights
        """
        order = ((window - 1) // 2) + 1
        nwts = 2 * order + 1
        w = np.zeros([nwts])
        n = nwts // 2
        w[n] = 2 * cutoff
        k = np.arange(1., n)
        sigma = np.sin(np.pi * k / n) * n / (np.pi * k)
        firstfactor = np.sin(2. * np.pi * cutoff * k) / (np.pi * k)
        w[n-1:0:-1] = firstfactor * sigma
        w[n+1:-1] = firstfactor * sigma
        w = w[1:-1]
        w = w / np.sum(w)
        return w

    l_filter = lweights(wi, co)

    if dim == 1:
        filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1))
        for i in range(len(field)-len(l_filter)+1):
            filtered[i] = np.sum(field[i:i+len(l_filter)] * l_filter)
    if dim == 2:
        filtered = np.zeros((np.ma.size(field, axis=0)-len(l_filter)+1,
                            np.ma.size(field, axis=1)))
        for i in range(len(field)-len(l_filter)+1):
            filtered[i] = np.sum((field[i:i+len(l_filter)].transpose() *
                                  l_filter).transpose(), axis=0)
    return filtered


def hovext(u, uc, ue, vT, vTc, vTe, uv, uvc, uve, l_filter):
    import matplotlib.pyplot as plt
    plt.figure()
    time = np.arange(1, 724/4, 0.25)
    timem = np.arange(6, 175, 0.25)
    timem = np.arange(6.5, 174, 0.25)
    up = ((np.max(u[:, 31, :32], axis=1) - np.max(uc[31, :32])) *
          np.sign(np.max(u[:, 31, :32], axis=1)) /
          abs(np.max(ue[31, :32]) - np.max(uc[31, :32])))
    uvp = ((np.min(uv[:, 18, :32], axis=1) - np.min(uvc[18, :32])) *
           np.sign(np.min(uv[:, 18, :32], axis=1)) /
           abs(np.min(uve[18, :32]) - np.min(uvc[18, :32])))
    vTp = ((np.min(vT[:, 31, :32], axis=1) - np.min(vTc[31, :32])) *
           np.sign(np.min(vT[:, 31, :32], axis=1)) /
           abs(np.min(vTe[31, :32]) - np.min(vTc[31, :32])))
    upm = np.zeros((676))
    uvpm = np.zeros((676))
    vTpm = np.zeros((676))
    for i in range(676):
        upm[i] = np.mean(up[i:i+44])
        uvpm[i] = np.mean(uvp[i:i+44])
        vTpm[i] = np.mean(vTp[i:i+44])
    upm = np.zeros((670))
    uvpm = np.zeros((670))
    vTpm = np.zeros((670))
    for i in range(670):
        upm[i] = np.sum(up[i:i+51] * l_filter)
        uvpm[i] = np.sum(uvp[i:i+51] * l_filter)
        vTpm[i] = np.sum(vTp[i:i+51] * l_filter)
    plt.plot(time, up, color='k', alpha=.5)
    plt.plot(time, uvp, color='r', alpha=.5)
    plt.plot(time, vTp, color='b', alpha=.5)
    plt.plot(timem, upm, color='k', label='u')
    plt.plot(timem, uvpm, color='r', label="u'v'")
    plt.plot(timem, vTpm, color='b', label="v'T'")
    plt.axhline(color='k')
    plt.axhline(-1, color='k')
    plt.xlabel('Time (days)')
    plt.xlim([0, 180])
    plt.legend()


def hovmollersu1(field, fieldcontrol, lat, heatind, var='u'):
    import matplotlib.pyplot as plt
    plt.figure()
    ctrs = {}
    ctrs['u'] = np.array([-3.2, -2.4, -1.6, -.8, .8, 1.6, 2.4, 3.2])
    ctrs['vT'] = np.array([-3.2, -2.4, -1.6, -.8, .8, 1.6, 2.4, 3.2]) * 1.875
    ctrs['uv'] = np.array([-3.2, -2.4, -1.6, -.8, .8, 1.6, 2.4, 3.2]) * 5
    ctrs['duv'] = np.array([-1, -0.75, -.5, -.25, .25, .5, .75, 1]) * 1.6
    ctrs['dvT'] = np.array([-1, -0.75, -.5, -.25, .25, .5, .75, 1]) * .6
    time = np.arange(1, 724/4, 0.25)
    timem = np.arange(6.5, 174, 0.25)
    mesht, meshlat = np.meshgrid(timem, lat)
    filt_res = lanczos(field-fieldcontrol, dim=2).transpose()
    # filt_res /= np.max(abs(filt_res))
    cs = plt.contour(mesht, meshlat, filt_res, ctrs[var], colors='k')
    if var == 'duv':
        plt.clabel(cs, fmt='%1.1f', inline_spacing=-2)
        plt.axhline(lat[np.argmin(fieldcontrol)], color='b', ls='--', lw=2,
                    dashes=[7, 3])
        plt.axhline(lat[np.argmax(fieldcontrol)], color='b', ls='--', lw=2,
                    dashes=[7, 3])
    elif var == 'dvT':
        plt.axhline(lat[np.argmin(fieldcontrol)], color='b', ls='--', lw=2,
                    dashes=[7, 3])
        plt.axhline(lat[np.argmax(fieldcontrol)], color='b', ls='--', lw=2,
                    dashes=[7, 3])
        plt.clabel(cs, fmt='%1.2f', inline_spacing=-2)
    else:
        plt.axhline(lat[np.argmax(abs(fieldcontrol))], color='b', ls='--',
                    lw=2, dashes=[7, 3])
        plt.clabel(cs, fmt='%1.1f', inline_spacing=-2)
    plt.axhline(lat[heatind], color='r', ls='--', lw=2,
                dashes=[7, 3])
    plt.xlabel('Time (days)')
    plt.ylabel('Latitude ($^\circ$)')
    plt.xlim([0, 180])
    plt.yticks(np.arange(-75, 90, 15))
    plt.ylim([lat[0], lat[-1]])
    plt.title(var)


def hovmollersu3(field, fieldcontrol, lat, heatind, var='u'):
    import matplotlib.pyplot as plt
    plt.figure()
    ctrs = {}
    ctrs['u'] = np.array([-3.2, -2.4, -1.6, -.8, .8, 1.6, 2.4, 3.2]) * 2
    ctrs['vT'] = np.array([-3.2, -2.4, -1.6, -.8, .8, 1.6, 2.4, 3.2]) * 1
    ctrs['uv'] = np.array([-3.2, -2.4, -1.6, -.8, .8, 1.6, 2.4, 3.2]) * 5
    ctrs['duv'] = np.array([-1, -0.75, -.5, -.25, .25, .5, .75, 1]) * 1.6
    ctrs['dvT'] = np.array([-1, -0.75, -.5, -.25, .25, .5, .75, 1]) * .4
    time = np.arange(1, 724/4, 0.25)
    mesht, meshlat = np.meshgrid(time, lat)
    cs = plt.contour(mesht, meshlat, (field-fieldcontrol).transpose(),
                     ctrs[var], colors='k')
    if var == 'duv':
        plt.clabel(cs, fmt='%1.1f', inline_spacing=-2)
        plt.axhline(lat[np.argmin(fieldcontrol[32:])+32], color='b', ls='--', lw=2,
                    dashes=[7, 3])
        plt.axhline(lat[np.argmax(fieldcontrol[32:])+32], color='b', ls='--', lw=2,
                    dashes=[7, 3])
    elif var == 'dvT':
        plt.axhline(lat[np.argmin(fieldcontrol[32:])+32], color='b', ls='--', lw=2,
                    dashes=[7, 3])
        plt.axhline(lat[np.argmax(fieldcontrol[32:])+32], color='b', ls='--', lw=2,
                    dashes=[7, 3])
        plt.clabel(cs, fmt='%1.2f', inline_spacing=-2)
    else:
        plt.axhline(lat[np.argmax(abs(fieldcontrol[32:]))+32], color='b', ls='--',
                    lw=2, dashes=[7, 3])
        plt.clabel(cs, fmt='%1.1f', inline_spacing=-2)
    plt.axhline(lat[heatind], color='r', ls='--', lw=2,
                dashes=[7, 3])
    plt.xlabel('Time (days)')
    plt.ylabel('Latitude ($^\circ$)')
    plt.xlim([0, 180])
    plt.yticks(np.arange(-75, 90, 15))
    plt.ylim([lat[0], lat[-1]])
    plt.title(var)
    print(np.min(field-fieldcontrol).transpose())


def hovmollersubw(u, uc, ue, vT, vTc, vTe, uv, uvc, uve, lat, pfull, heatind,
                  maxlat=90):
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines

    l1 = 18
    l2 = 24
    l3 = 31
    duv = diverg(uv, lat, pfull)
    duvc = diverg(np.expand_dims(uvc, axis=0), lat, pfull)[0]
    duve = diverg(np.expand_dims(uve, axis=0), lat, pfull)[0]

    duvp1 = ((np.min(duv[:, l1, :32], axis=1) - np.min(duvc[l1, :32])) *
             np.sign(np.min(duv[:, l1, :32], axis=1)) /
             abs(np.min(duve[l1, :32]) - np.min(duvc[l1, :32])))
    up3 = ((np.max(u[:, l3, :32], axis=1) - np.max(uc[l3, :32])) *
           np.sign(np.max(u[:, l3, :32], axis=1)) /
           abs(np.max(ue[l3, :32]) - np.max(uc[l3, :32])))
    vTp3 = ((np.min(vT[:, l3, :32], axis=1) - np.min(vTc[l3, :32])) *
            np.sign(np.min(vT[:, l3, :32], axis=1)) /
            abs(np.min(vTe[l3, :32]) - np.min(vTc[l3, :32])))

    up3 = lanczos(up3)
    duvp1 = lanczos(duvp1)
    vTp3 = lanczos(vTp3)
    u2 = lanczos(u[:, l2]-uc[l2], dim=2).transpose()
    u2 /= np.max(abs(u2))
    u3 = lanczos(u[:, l3]-uc[l3], dim=2).transpose()
    u3 /= np.max(abs(u3))
    duv1 = lanczos(duv[:, l1]-duvc[l1], dim=2).transpose()
    duv1 /= np.max(abs(duv1))
    duv2 = lanczos(duv[:, l2]-duvc[l2], dim=2).transpose()
    duv2 /= np.max(abs(duv2))
    vT2 = lanczos(vT[:, l2]-vTc[l2], dim=2).transpose()
    vT2 /= np.max(abs(vT2))
    vT3 = lanczos(vT[:, l3]-vTc[l3], dim=2).transpose()
    vT3 /= np.max(abs(vT3))

    timem = np.arange(6.5, 174, 0.25)
    ctrs = [-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1]
    mesht, meshlat = np.meshgrid(timem, lat)
    red_line = mlines.Line2D([], [], color='r', linestyle='-',
                             linewidth=2, alpha=0.75)
    blue_line = mlines.Line2D([], [], color='b', linestyle='--',
                              linewidth=2, dashes=(2.5, 2.5), alpha=0.75)
    black_line = mlines.Line2D([], [], color='k', linestyle='-',
                               linewidth=2)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    ax1 = fig.add_subplot(2, 2, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.plot(timem, vTp3, color='b', label="$v'T'_{850}$")
    ax1.plot(timem, duvp1, color='r', label=r"$\bigtriangledown\cdot u'v'_{250}$")
    ax1.plot(timem, up3, color='k', label='$u_{850}$')
    ax1.axhline(color='k', lw=2)
    ax1.axhline(-1, color='k', lw=2, ls='--')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('$\Delta\phi$')
    ax1.set_xlim([0, 173])
    ax1.set_yticks(np.arange(-4, 16, .5))
    ax1.set_ylim([-1.5, .5])
    ax1.legend()
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)

    ax2 = fig.add_subplot(2, 2, 2)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.contourf(mesht, meshlat, vT2, ctrs, cmap='bwr')
    cs = ax2.contour(mesht, meshlat, vT3, ctrs, colors='k')
    ax2.axhline(lat[np.argmin(vTc[l3])], color='g', ls='-', lw=2)
    ax2.axhline(lat[np.argmax(vTc[l3, 32:])+32], color='g', ls='-', lw=2)
    ax2.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax2.axhline(lat[heatind], color='blue', ls='-', lw=2)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Latitude ($^\circ$)')
    ax2.set_xlim([0, 173])
    ax2.set_yticks(np.arange(-90, 115, 15))
    ax2.set_ylim([-90, maxlat])
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)
    ax2.legend(handles=[(red_line, blue_line), black_line],
               labels=["$v'T'_{500}$", "$v'T'_{850}$"])

    ax3 = fig.add_subplot(2, 2, 3)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax3.spines[axis].set_linewidth(2)
    ax3.contourf(mesht, meshlat, duv2, ctrs, cmap='bwr')
    cs = ax3.contour(mesht, meshlat, duv1, ctrs, colors='k')
    ax3.axhline(lat[np.argmin(duvc[l1])], color='g', ls='-', lw=2)
    ax3.axhline(lat[np.argmin(duvc[l1, 32:])+32], color='g', ls='-', lw=2)
    ax3.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax3.axhline(lat[heatind], color='blue', ls='-', lw=2)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Latitude ($^\circ$)')
    ax3.set_xlim([0, 173])
    ax3.set_yticks(np.arange(-90, 115, 15))
    ax3.set_ylim([-90, maxlat])
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)
    ax3.legend(handles=[(red_line, blue_line), black_line],
               labels=[r"$\bigtriangledown\cdot u'v'_{500}$",
                       r"$\bigtriangledown\cdot u'v'_{250}$"])

    ax4 = fig.add_subplot(2, 2, 4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax4.spines[axis].set_linewidth(2)
    ax4.contourf(mesht, meshlat, u2, ctrs, cmap='bwr')
    cs = ax4.contour(mesht, meshlat, u3, ctrs, colors='k')
    ax4.axhline(lat[np.argmax(uc[l3])], color='g', ls='-', lw=2)
    ax4.axhline(lat[np.argmax(uc[l3, 32:])+32], color='g', ls='-', lw=2)
    ax4.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax4.axhline(lat[heatind], color='blue', ls='-', lw=2)
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Latitude ($^\circ$)')
    ax4.set_xlim([0, 173])
    ax4.set_yticks(np.arange(-90, 115, 15))
    ax4.set_ylim([-90, maxlat])
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)
    ax4.legend(handles=[(red_line, blue_line), black_line],
               labels=['$u_{500}$', '$u_{850}$'])

    plt.subplots_adjust(hspace=0.15, wspace=0.15*3/4, top=.95, bottom=0.05,
                        left=0.05, right=.99)


def hovmollersubs(u, uc, ue, vT, vTc, vTe, uv, uvc, uve, lat, pfull, heatind):
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    latmin = 0
    l1 = 18
    l2 = 24
    l3 = 31
    duv = diverg(uv, lat, pfull)
    duvc = diverg(np.expand_dims(uvc, axis=0), lat, pfull)[0]
    duve = diverg(np.expand_dims(uve, axis=0), lat, pfull)[0]

    duvp1 = ((np.min(duv[:, l1, 32:], axis=1) - np.min(duvc[l1, 32:])) *
             np.sign(np.min(duv[:, l1, 32:], axis=1)) /
             abs(np.min(duve[l1, 32:]) - np.min(duvc[l1, 32:])))
    up3 = ((np.max(u[:, l3, 32:], axis=1) - np.max(uc[l3, 32:])) *
           np.sign(np.max(u[:, l3, 32:], axis=1)) /
           abs(np.max(ue[l3, 32:]) - np.max(uc[l3, 32:])))
    vTp3 = ((np.max(vT[:, l3, 32:], axis=1) - np.max(vTc[l3, 32:])) *
            np.sign(np.max(vT[:, l3, 32:], axis=1)) /
            abs(np.max(vTe[l3, 32:]) - np.max(vTc[l3, 32:])))

    up3 = lanczos(up3)
    duvp1 = lanczos(duvp1)
    vTp3 = lanczos(vTp3)
    u2 = lanczos(u[:, l2]-uc[l2], dim=2).transpose()
    u2 /= np.max(abs(u2))
    u3 = lanczos(u[:, l3]-uc[l3], dim=2).transpose()
    u3 /= np.max(abs(u3))
    duv1 = lanczos(duv[:, l1]-duvc[l1], dim=2).transpose()
    duv1 /= np.max(abs(duv1))
    duv2 = lanczos(duv[:, l2]-duvc[l2], dim=2).transpose()
    duv2 /= np.max(abs(duv2))
    vT2 = lanczos(vT[:, l2]-vTc[l2], dim=2).transpose()
    vT2 /= np.max(abs(vT2))
    vT3 = lanczos(vT[:, l3]-vTc[l3], dim=2).transpose()
    vT3 /= np.max(abs(vT3))

    timem = np.arange(6.5, 174, 0.25)
    ctrs = [-1, -.8, -.6, -.4, -.2, .2, .4, .6, .8, 1]
    mesht, meshlat = np.meshgrid(timem, lat)
    red_line = mlines.Line2D([], [], color='r', linestyle='-',
                             linewidth=2, alpha=0.75)
    blue_line = mlines.Line2D([], [], color='b', linestyle='--',
                              linewidth=2, dashes=(2.5, 2.5), alpha=0.75)
    black_line = mlines.Line2D([], [], color='k', linestyle='-',
                               linewidth=2)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)
    ax1 = fig.add_subplot(2, 2, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.plot(timem, vTp3, color='b', label="$v'T'_{850}$")
    ax1.plot(timem, duvp1, color='r', label=r"$\bigtriangledown\cdot u'v'_{250}$")
    ax1.plot(timem, up3, color='k', label='$u_{850}$')
    ax1.axhline(color='k', lw=2)
    ax1.axhline(1, color='k', lw=2, ls='--')
    ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('$\Delta\phi$')
    ax1.set_xlim([0, 173])
    ax1.set_yticks(np.arange(-2, 16, 2))
    ax1.set_ylim([0, 12])
    ax1.legend()
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)

    ax2 = fig.add_subplot(2, 2, 2)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.contourf(mesht, meshlat, vT2, ctrs, cmap='bwr')
    cs = ax2.contour(mesht, meshlat, vT3, ctrs, colors='k')
    ax2.axhline(lat[np.argmax(vTc[l3, 32:])+32], color='g', ls='-', lw=2)
    ax2.axhline(lat[np.argmin(vTc[l3])], color='g', ls='-', lw=2)
    ax2.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax2.axhline(lat[heatind], color='blue', ls='-', lw=2)
    ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Latitude ($^\circ$)')
    ax2.set_xlim([0, 173])
    ax2.set_yticks(np.arange(-90, 115, 15))
    ax2.set_ylim([latmin, 90])
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)
    ax2.legend(handles=[(red_line, blue_line), black_line],
               labels=["$v'T'_{500}$", "$v'T'_{850}$"])

    ax3 = fig.add_subplot(2, 2, 3)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax3.spines[axis].set_linewidth(2)
    ax3.contourf(mesht, meshlat, duv2, ctrs, cmap='bwr')
    cs = ax3.contour(mesht, meshlat, duv1, ctrs, colors='k')
    ax3.axhline(lat[np.argmin(duvc[l1, 32:])+32], color='g', ls='-', lw=2)
    ax3.axhline(lat[np.argmin(duvc[l1])], color='g', ls='-', lw=2)
    ax3.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax3.axhline(lat[heatind], color='blue', ls='-', lw=2)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Latitude ($^\circ$)')
    ax3.set_xlim([0, 173])
    ax3.set_yticks(np.arange(-90, 115, 15))
    ax3.set_ylim([latmin, 90])
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)
    ax3.legend(handles=[(red_line, blue_line), black_line],
               labels=[r"$\bigtriangledown\cdot u'v'_{500}$",
                       r"$\bigtriangledown\cdot u'v'_{250}$"])

    ax4 = fig.add_subplot(2, 2, 4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax4.spines[axis].set_linewidth(2)
    ax4.contourf(mesht, meshlat, u2, ctrs, cmap='bwr')
    cs = ax4.contour(mesht, meshlat, u3, ctrs, colors='k')
    ax4.axhline(lat[np.argmax(uc[l3, 32:])+32], color='g', ls='-', lw=2)
    ax4.axhline(lat[np.argmax(uc[l3])], color='g', ls='-', lw=2)
    ax4.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax4.axhline(lat[heatind], color='blue', ls='-', lw=2)
    ax4.set_xlabel('Time (days)')
    ax4.set_ylabel('Latitude ($^\circ$)')
    ax4.set_xlim([0, 173])
    ax4.set_yticks(np.arange(-90, 115, 15))
    ax4.set_ylim([latmin, 90])
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.075, y=0.95)
    ax4.legend(handles=[(red_line, blue_line), black_line],
               labels=['$u_{500}$', '$u_{850}$'])

    plt.subplots_adjust(hspace=0.15, wspace=0.15*3/4, top=.95, bottom=0.05,
                        left=0.05, right=.99)


def animatesu(data, data2, potempcontrol, temp, z, sigma, lat):
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
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    my_cmap = plt.cm.RdYlBu_r(np.arange(256))
    my_cmap[120:136, :] = 0.99
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)

    #  U and pottemp
    ctrs1 = np.arange(-10, 0, 1)
    ctrs2 = np.arange(1, 11, 1)
    ctrs = np.concatenate((ctrs1, ctrs2))
    ctrs3 = [270, 280, 290, 300, 310, 320, 330, 340, 350, 400, 500]

    # vT and uv
    ctrs3 = np.array([-20, -15, -10, -5, 5, 10, 15, 20, 25, 30])*5e-7
    ctrs = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
    fig = plt.figure()

    def updatefig(i):
        fig.clear()
        lapse0 = lapserate(temp[0, :, :], z[0, :, :], sigma, lat)
        lapse0[22:, :] = 99
        lapse = lapserate(temp[i, :, :], z[i, :, :], sigma, lat)
        lapse[22:, :] = 99
        plot = plt.contourf(meshlat, meshsigma, data[i, :, :], ctrs,
                            cmap=newcmap, vmin=-10, vmax=10, extend='both')
        plt.contour(meshlat, meshsigma, potempcontrol,
                    ctrs3, colors='gray', linewidths=1.5)
        plot1 = plt.contour(meshlat, meshsigma, data2[i, :, :], ctrs3,
                            colors='k', linewidths=1.5)
        plt.contour(meshlat[1:, :], meshsigma[1:, :], lapse0, [2],
                    linewidths=2, colors='dodgerblue')
        plt.contour(meshlat[1:, :], meshsigma[1:, :], lapse, [2], linewidths=2)
        # plt.scatter(lat[20], sigma[22], marker='o', color='g',
        # plt.scatter(lat[8], sigma[33], marker='o', color='g',
        # plt.scatter(lat[24], sigma[22], marker='o', color='g',
        # plt.scatter(lat[49], sigma[33], marker='o', color='g',
        plt.scatter(lat[41], sigma[16], marker='o', color='g',
                    linewidth=2, s=50)
        plt.clabel(plot1, inline=True, inline_spacing=-3, fontsize=10,
                   fmt='%.0f')
        plt.gca().invert_yaxis()
        plt.yscale('linear')
        plt.xlim([-90, 90])
        plt.xticks(np.arange(-90, 105, 15))
        plt.title('Day ' + str((i+1)/4), y=1.08, fontsize=30)
        plt.colorbar(plot, orientation='horizontal',
                     shrink=0.5, spacing='proportional',
                     label="Change in vT [K m/s] (Contours show change in u'v')")
        plt.draw()

    anim = animation.FuncAnimation(fig, updatefig, np.ma.size(data, axis=0))
    #  anim.save('test.mp4',bitrate=10000)
    return anim


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
