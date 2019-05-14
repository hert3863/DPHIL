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
    for a in range(len(pfull)):
        pottemp[a] = temp[a] * (1000 / pfull) ** k
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


def epcalc(uv, vT, T, lat, pfull):
    if uv.ndim == 2:
        uv = np.expand_dims(uv, axis=0)
        vT = np.expand_dims(vT, axis=0)
        T = np.expand_dims(T, axis=0)
    k = 0.286
    w = 7.2921e-5
    a = 6.371e6
    meshlat, meshp = np.meshgrid(lat, pfull)
    meshlatc = np.cos(np.pi*meshlat/180)
    meshlats = np.sin(np.pi*meshlat/180)
    f = 2 * w * meshlats
    vTheta = vT * (1000 / meshp) ** k
    theta = T * (1000 / meshp) ** k
    stat_stab = np.zeros(np.shape(theta))
    lpfull = np.log(pfull*100)
    stat_stab[:, 0, :] = (0.01/pfull[0]) * ((theta[:, 1, :] - theta[:, 0, :]) /
                                            (lpfull[1]-lpfull[0]))
    stat_stab[:, -1, :] = (0.01/pfull[-1]) * ((theta[:, -1, :] -
                                               theta[:, -2, :]) /
                                              (lpfull[-1]-lpfull[-2]))
    for i in range(len(pfull)-2):
        stat_stab[:, i+1, :] = (0.01/pfull[i+1]) * ((theta[:, i+2, :] -
                                                    theta[:, i, :]) /
                                                    (lpfull[i+2]-lpfull[i]))
    ep_p = f * a * vTheta / stat_stab
    ep_y = -a * meshlatc * uv

    dep_p = np.zeros(np.shape(ep_p))
    dep_p[:, 0, :] = (0.01/pfull[0]) * ((ep_p[:, 1, :] - ep_p[:, 0, :]) /
                                        (lpfull[1]-lpfull[0]))
    dep_p[:, -1, :] = (0.01/pfull[-1]) * ((ep_p[:, -1, :] - ep_p[:, -2, :]) /
                                          (lpfull[-1]-lpfull[-2]))
    for i in range(len(pfull)-2):
        dep_p[:, i+1, :] = (0.01/pfull[i+1]) * ((ep_p[:, i+2, :] -
                                                ep_p[:, i, :]) /
                                                (lpfull[i+2]-lpfull[i]))

    dep_y = np.zeros(np.shape(ep_y))
    aslat = a * np.sin(lat*np.pi/180)
    dep_y[:, :, 0] = ((ep_y[:, :, 1] - ep_y[:, :, 0]) /
                      (aslat[1]-aslat[0]))
    dep_y[:, :, -1] = ((ep_y[:, :, -1] - ep_y[:, :, -2]) /
                       (aslat[-1]-aslat[-2]))
    for i in range(len(lat)-2):
        dep_y[:, :, i+1] = ((ep_y[:, :, i+2] - ep_y[:, :, i]) /
                            (aslat[i+2]-aslat[i]))

    divF = (dep_y + dep_p) / (a*meshlatc/86400)

    ep_y = ep_y * np.sqrt(1000/meshp) / (np.pi*a)
    ep_p = ep_p * meshlatc * np.sqrt(1000/meshp) / 1e5
    return ep_y, ep_p, divF


def epplot(uv, vT, temp, uvc, vTc, tempc, uve, vTe, tempe, lat, pfull, heat,
           minlat=-90, maxlat=90, hid=0, dc=8, d=8, mc=1, m=.5, ld=1):
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    ep_y, ep_p, divF = epcalc(uv, vT, temp, lat, pfull)
    ep_yc, ep_pc, divFc = epcalc(uvc, vTc, tempc, lat, pfull)
    ep_ye, ep_pe, divFe = epcalc(uve[hid], vTe[hid], tempe[hid], lat, pfull)
    duv = diverg(uv, lat, pfull) * 86400
    duvc = diverg(uvc, lat, pfull) * 86400
    duve = diverg(uve[hid], lat, pfull) * 86400
    meshlat, meshp = np.meshgrid(lat, pfull)
    ctrs = np.linspace(-1, 1, 17)
    ctrs2 = np.array([-4, -3, -2, -1, 1, 2, 3, 4])
    sigmax = np.argmax(heat[hid, :, :], axis=0)[0]
    if sigmax == 33:
        sigmax = 31
    latmax = np.argmax(heat[hid, :, :], axis=1)[0]
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)

    fig = plt.figure(facecolor='w', edgecolor='k', linewidth=2)

    ax1 = fig.add_subplot(2, 2, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.contourf(meshlat, meshp, divFc[0], ctrs*dc,
                 cmap=newcmap, extend='both')
    ax1.quiver(meshlat[1:32:2, ::ld], meshp[1:32:2, ::ld], ep_yc[0, 1:32:2, ::ld], -ep_pc[0, 1:32:2, ::ld],
               scale=1000, angles='uv', width=0.0075, minlength=0.5)
    ax1.contour(meshlat, meshp, duvc[0], ctrs2*mc,
                colors=('r', 'r', 'r', 'r', 'k', 'k', 'k', 'k'))
    ax1.scatter(lat[latmax], pfull[sigmax], marker='o',
                color='deeppink', linewidth=2, s=50)
    # ax1.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax1.set_ylabel('Pressure (hPa)', fontsize=16)
    ax1.set_ylim([1000, 50])
    ax1.set_xticks(np.arange(-90, 115, 15))
    ax1.set_xlim([minlat, maxlat])
    ax1.tick_params(axis='both', labelsize=14)
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.045, y=0.95)

    ax2 = fig.add_subplot(2, 2, 2)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.contourf(meshlat, meshp, np.mean(divF[12:32], axis=0)-divFc[0], ctrs*d,
                 cmap=newcmap, extend='both')
    ax2.quiver(meshlat[1:32:2, ::ld], meshp[1:32:2, ::ld],
               np.mean(ep_y[12:32, 1:32:2, ::ld], axis=0)-ep_yc[0, 1:32:2, ::ld],
               np.mean(-ep_p[12:32, 1:32:2, ::ld], axis=0)+ep_pc[0, 1:32:2, ::ld],
               scale=500, angles='uv', width=0.0075, minlength=0.5)
    ax2.contour(meshlat, meshp, np.mean(duv[12:32], axis=0)-duvc[0],
                ctrs2*m, colors=('r', 'r', 'r', 'r', 'k', 'k', 'k', 'k'))
    ax2.scatter(lat[latmax], pfull[sigmax], marker='o',
                color='deeppink', linewidth=2, s=50)
    ax2.set_ylim([850, 50])
    ax2.set_xticks(np.arange(-90, 115, 15))
    ax2.set_xlim([minlat, maxlat])
    ax2.tick_params(axis='both', labelsize=14)
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.045, y=0.95)

    di1 = 108
    di2 = 128

    ax3 = fig.add_subplot(2, 2, 3)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax3.spines[axis].set_linewidth(2)
    ax3.contourf(meshlat, meshp, np.mean(divF[di1:di2], axis=0)-np.mean(divF[12:32], axis=0), ctrs*d,
                 cmap=newcmap, extend='both')
    ax3.quiver(meshlat[1:32:2, ::ld], meshp[1:32:2, ::ld],
               np.mean(ep_y[di1:di2, 1:32:2, ::ld], axis=0)-np.mean(ep_y[12:32, 1:32:2, ::ld], axis=0),
               np.mean(-ep_p[di1:di2, 1:32:2, ::ld], axis=0)-np.mean(-ep_p[12:32, 1:32:2, ::ld], axis=0),
               scale=500, angles='uv', width=0.0075, minlength=0.5)
    ax3.contour(meshlat, meshp, np.mean(duv[di1:di2], axis=0)-np.mean(duv[12:32], axis=0),
                ctrs2*m, colors=('r', 'r', 'r', 'r', 'k', 'k', 'k', 'k'))
    ax3.scatter(lat[latmax], pfull[sigmax], marker='o',
                color='deeppink', linewidth=2, s=50)
    ax3.set_ylabel('Pressure (hPa)', fontsize=16)
    ax3.set_xlabel('Latitude ($^\circ$)', fontsize=16)
    ax3.set_ylim([850, 50])
    ax3.set_xticks(np.arange(-90, 115, 15))
    ax3.set_xlim([minlat, maxlat])
    ax3.tick_params(axis='both', labelsize=14)
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.045, y=0.95)

    ax4 = fig.add_subplot(2, 2, 4)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax4.spines[axis].set_linewidth(2)
    plot1 = ax4.contourf(meshlat, meshp, divFe[0]-np.mean(divF[di1:di2], axis=0), ctrs*d,
                         cmap=newcmap, extend='both')
    ax4.quiver(meshlat[1:32:2, ::ld], meshp[1:32:2, ::ld], ep_ye[0, 1:32:2, ::ld]-np.mean(ep_y[di1:di2, 1:32:2, ::ld], axis=0),
               -ep_pe[0, 1:32:2, ::ld]-np.mean(-ep_p[di1:di2, 1:32:2, ::ld], axis=0),
               scale=500, angles='uv', width=0.0075, minlength=0.5)
    ax4.contour(meshlat, meshp, duve[0]-np.mean(duv[di1:di2], axis=0), ctrs2*m,
                colors=('r', 'r', 'r', 'r', 'k', 'k', 'k', 'k'))
    ax4.scatter(lat[latmax], pfull[sigmax], marker='o',
                color='deeppink', linewidth=2, s=50)
    ax4.set_xlabel('Latitude ($^\circ$)', fontsize=16)
    ax4.set_ylim([850, 50])
    ax4.set_xticks(np.arange(-90, 115, 15))
    ax4.set_xlim([minlat, maxlat])
    ax4.tick_params(axis='both', labelsize=14)
    ax4.set_title('d', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.045, y=0.95)

    cbar_ax = fig.add_axes([0.35, 0.065, 0.35, 0.01])
    b = fig.colorbar(plot1, cax=cbar_ax, spacing='proportional',
                     orientation='horizontal', extend='max', format='%.1f')
    b.set_label(label='EP flux divergence (ms$^{-1}$ day$^{-1}$)', size=12,
                fontsize=16, fontname='Arial', labelpad=-2.5)
    cl = plt.getp(cbar_ax, 'xmajorticklabels')
    plt.setp(cl, fontname='Arial', fontsize=16)
    for label in cbar_ax.xaxis.get_ticklabels()[1::2]:
        label.set_visible(False)

    ax1.annotate('Ctrl', xy=(0.5, 1), xytext=(0, 6), fontsize=16,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax2.annotate('Day 3-7 minus ctrl', xy=(0.5, 1), xytext=(0, 6), fontsize=16,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax3.annotate('Day 28-32 minus day 3-7', xy=(0.5, 1), xytext=(0, 6), fontsize=16,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    ax4.annotate('Eqm minus day 28-32', xy=(0.5, 1), xytext=(0, 6), fontsize=16,
                 fontname='Arial',
                 xycoords='axes fraction', textcoords='offset points',
                 ha='center', va='baseline')
    plt.subplots_adjust(wspace=0.11, hspace=0.18, top=.97, bottom=0.12,
                        left=0.075, right=.975)


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
        u = np.transpose(np.transpose(u) - np.transpose(np.mean(u, axis=3)))
        v = np.transpose(np.transpose(v) - np.transpose(np.mean(v, axis=3)))
        T = np.transpose(np.transpose(T) - np.transpose(np.mean(T, axis=3)))
        #uv1 = (u - np.mean(u, axis=0)) * (v - np.mean(v, axis=0))
        #vT1 = (T - np.mean(T, axis=0)) * (v - np.mean(v, axis=0))
        uv1 = u * v
        vT1 = T * v
        uv1 = np.mean(uv1, axis=3)
        vT1 = np.mean(vT1, axis=3)
        uv1 = np.mean(uv1, axis=0)
        vT1 = np.mean(vT1, axis=0)
        uv[i, :] = uv1
        vT[i, :] = vT1
        print(str(i))
    return uv, vT


def diverg(field, lat, pfull):
    if field.ndim == 2:
        field = np.expand_dims(field, axis=0)
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


def hovmollersub(u, uc, vT, vTc, uv, uvc, lat, pfull, heatind, minlat=90,
                 maxlat=90):
    import matplotlib.pyplot as plt
    import matplotlib.lines as mlines
    l1 = 18
    l2 = 24
    l3 = 31
    duv = diverg(uv, lat, pfull)
    duvc = diverg(np.expand_dims(uvc, axis=0), lat, pfull)[0]

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

    ax1 = fig.add_subplot(3, 1, 1)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax1.spines[axis].set_linewidth(2)
    ax1.contourf(mesht, meshlat, vT2, ctrs, cmap='bwr')
    cs = ax1.contour(mesht, meshlat, vT3, ctrs, colors='k')
    ax1.axhline(lat[np.argmin(vTc[l3])], color='g', ls='-', lw=2)
    ax1.axhline(lat[np.argmax(vTc[l3, 32:])+32], color='g', ls='-', lw=2)
    ax1.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax1.axhline(lat[heatind], color='deeppink', ls='-', lw=2)
    # ax1.set_xlabel('Time (days)')
    ax1.set_ylabel('Latitude ($^\circ$)')
    ax1.set_xlim([0, 173])
    ax1.set_yticks(np.arange(-90, 115, 15))
    ax1.set_ylim([minlat, maxlat])
    ax1.set_title('a', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.95)
    ax1.legend(handles=[(red_line, blue_line), black_line],
               labels=["$v'T'_{500}$", "$v'T'_{850}$"], loc=1)

    ax2 = fig.add_subplot(3, 1, 2)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax2.spines[axis].set_linewidth(2)
    ax2.contourf(mesht, meshlat, duv2, ctrs, cmap='bwr')
    cs = ax2.contour(mesht, meshlat, duv1, ctrs, colors='k')
    ax2.axhline(lat[np.argmin(duvc[l1, 32:])+32], color='g', ls='-', lw=2)
    ax2.axhline(lat[np.argmin(duvc[l1])], color='g', ls='-', lw=2)
    ax2.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax2.axhline(lat[heatind], color='deeppink', ls='-', lw=2)
    # ax2.set_xlabel('Time (days)')
    ax2.set_ylabel('Latitude ($^\circ$)')
    ax2.set_xlim([0, 173])
    ax2.set_yticks(np.arange(-90, 115, 15))
    ax2.set_ylim([minlat, maxlat])
    ax2.set_title('b', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.95)
    ax2.legend(handles=[(red_line, blue_line), black_line],
               labels=[r"$\bigtriangledown\cdot u'v'_{500}$",
                       r"$\bigtriangledown\cdot u'v'_{250}$"], loc=1)

    ax3 = fig.add_subplot(3, 1, 3)
    for axis in ['top', 'bottom', 'left', 'right']:
        ax3.spines[axis].set_linewidth(2)
    ax3.contourf(mesht, meshlat, u2, ctrs, cmap='bwr')
    cs = ax3.contour(mesht, meshlat, u3, ctrs, colors='k')
    ax3.axhline(lat[np.argmax(uc[l3, 32:])+32], color='g', ls='-', lw=2)
    ax3.axhline(lat[np.argmax(uc[l3])], color='g', ls='-', lw=2)
    ax3.clabel(cs, fmt='%1.1f', inline_spacing=-6)
    ax3.axhline(lat[heatind], color='deeppink', ls='-', lw=2)
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('Latitude ($^\circ$)')
    ax3.set_xlim([0, 173])
    ax3.set_yticks(np.arange(-90, 115, 15))
    ax3.set_ylim([minlat, maxlat])
    ax3.set_title('c', loc='left', fontname='Arial', fontsize=20,
                  fontweight='bold', x=-.1, y=0.95)
    ax3.legend(handles=[(red_line, blue_line), black_line],
               labels=['$u_{500}$', '$u_{850}$'], loc=1)

    plt.subplots_adjust(hspace=0.1, top=.99, bottom=0.05,
                        left=0.1, right=.99)


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


def postage(potemp, potempcontrol, vT, lat, pfull):
    import numpy as np
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    meshlat, meshp = np.meshgrid(lat, pfull)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("nj", my_cmap)
    ctrs = np.arange(-2, 2.2, .2) 
    ctrs3 = np.arange(270, 350, 10)

    fig, axs = plt.subplots(4, 5, facecolor='w', edgecolor='k', sharex=True,
                            sharey=True, linewidth=2)

    for i in range(20):
        ax = axs[int(np.floor(i/5)), np.remainder(i, 5)]
        for axis in ['top', 'bottom', 'left', 'right']:
            ax.spines[axis].set_linewidth(2)
        plot = ax.contourf(meshlat, meshp, vT[4*i, :, :]-vT[0], ctrs,
                           cmap=newcmap, vmin=np.min(ctrs), vmax=np.max(ctrs),
                           extend='both')
        ax.contour(meshlat, meshp, potempcontrol,
                   ctrs3, colors='green', linewidths=1.5)
        plot1 = ax.contour(meshlat, meshp, potemp[4*i],
                           ctrs3, colors='black', linewidths=1.5)
        ax.scatter(lat[49], pfull[33], marker='o', color='deeppink',
                   linewidth=2, s=50)
        ax.clabel(plot1, inline=True, inline_spacing=-7,
                  fontsize=10, fmt='%.0f')
        plt.yscale('linear')
        plt.xticks(np.arange(-90, 105, 15))
        plt.xlim([0, 90])
        plt.ylim([1000, 100])
    #plt.title('Day ' + str((i+1)/4), y=1.08, fontsize=30)
    #plt.colorbar(plot, orientation='horizontal',
            #     shrink=0.5, spacing='proportional',
           #      label="v'T' response (K ms$^{-1}$)")
    plt.subplots_adjust(hspace=0.1, top=.99, bottom=0.05,
                        left=0.1, right=.99)


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
