# -*- coding: utf-8 -*-
"""
Created on Thu Sep 14 13:35:13 2017

@author: bakerh
"""

import matplotlib.pyplot as plt
import scipy.stats as stats
import numpy as np


def scatter(d_means, dT, hemisphere='Winter', cell='Hadley'):
    if hemisphere == 'Winter':
        s = 0
    else:
        s = 1
    plt.figure(figsize=(9, 9))
    a_s = stats.linregress(np.ndarray.flatten(d_means[s,:,17:]),np.ndarray.flatten(dT[s,:,17:]))
    a_w = stats.linregress(np.ndarray.flatten(d_means[s,:,:17]),np.ndarray.flatten(dT[s,:,:17]))
    plt.scatter(d_means[s,:,:17],dT[s,:,:17],color='b')
    plt.scatter(d_means[s,:,17:],dT[s,:,17:],color='r')
    plt.plot(np.arange(np.min(d_means[s]), np.max(d_means[s]),0.01),np.arange(np.min(d_means[s]), np.max(d_means[s]),0.01)*a_w[0]+a_w[1],color='b')
    plt.plot(np.arange(np.min(d_means[s]), np.max(d_means[s]),0.01),np.arange(np.min(d_means[s]), np.max(d_means[s]),0.01)*a_s[0]+a_s[1],color='r')
    plt.ylim([np.min(dT[s]), np.max(dT[s])])
    plt.xlim([np.min(d_means[s]), np.max(d_means[s])])
    plt.text(np.min(d_means[s])+0.01*abs(np.ptp(d_means[s])),
                    np.max(dT[s])-0.05*abs(np.ptp(dT[s])),'Summer: $\mathregular{r^2}=$ %.2f' % a_s[2]**2)
    plt.text(np.min(d_means[s])+0.01*abs(np.ptp(d_means[s])),
                    np.max(dT[s])-0.08*abs(np.ptp(dT[s])),'Winter: $\mathregular{r^2}=$ %.2f' % a_w[2]**2)
    plt.xlabel('Change in Hadley Cell extent (deg)')
    plt.ylabel('Change in heat flux across Hadley Cell')
    plt.title(hemisphere)


def scatter_sub(d_had, d_jet, eq_l):
    fig, axs = plt.subplots(2, 2, figsize=(11, 11),
                            facecolor='w', edgecolor='k', linewidth=2)
    for s in range(2):
        for j in range(2):
            a_s = stats.linregress(np.ndarray.flatten(d_had[1, s,:,eq_l:]), np.ndarray.flatten(d_jet[s, j,:,eq_l:]))
            a_w = stats.linregress(np.ndarray.flatten(d_had[1, s,:,:eq_l]), np.ndarray.flatten(d_jet[s, j,:,:eq_l]))
            axs[s, 1-j].scatter(d_had[1, s,:,:eq_l],d_jet[s, j,:,:eq_l],color='b')
            axs[s, 1-j].scatter(d_had[1, s,:,eq_l:],d_jet[s, j,:,eq_l:],color='r')
            axs[s, 1-j].plot(np.arange(np.min(d_had[1, s]), np.max(d_had[1, s]),0.01),np.arange(np.min(d_had[1, s]), np.max(d_had[1, s]),0.01)*a_w[0]+a_w[1],color='b')
            axs[s, 1-j].plot(np.arange(np.min(d_had[1, s]), np.max(d_had[1, s]),0.01),np.arange(np.min(d_had[1, s]), np.max(d_had[1, s]),0.01)*a_s[0]+a_s[1],color='r')
            if j == 1:
                axs[s, 1-j].plot(np.linspace(-6, 6, 100), np.linspace(-6, 6, 100), color='k', linestyle='--')
            axs[s, 1-j].set_ylim([np.min(d_jet[s, j]), np.max(d_jet[s, j])])
            axs[s, 1-j].set_xlim([np.min(d_had[1, s]), np.max(d_had[1, s])])
            axs[s, 1-j].text(np.max(d_had[1, s])-.35*abs(np.ptp(d_had[1, s])),
                            np.min(d_jet[s, j])+.05*abs(np.ptp(d_jet[s, j])),'Summer: $\mathregular{r^2}=$ %.2f' % a_s[2]**2)
            axs[s, 1-j].text(np.max(d_had[1, s])-0.32*abs(np.ptp(d_had[1, s])),
                            np.min(d_jet[s, j])+0.01*abs(np.ptp(d_jet[s, j])),'Winter: $\mathregular{r^2}=$ %.2f' % a_w[2]**2)
            axs[s, 1-j].set_xlabel('Poleward shift of Hadley cell extent (deg)')
            if j == 0:
                axs[s, 1-j].set_ylabel('Change in jet speed (ms$\mathregular{^{-1}}$)')
                if s == 0:
                    axs[s, 1-j].set_title('Magnitude')
            if j == 1:
                axs[s, 1-j].set_ylabel('Poleward shift of jet (deg)')
                if s == 0:
                    axs[s, 1-j].set_title('Latitude')
            # axs[s, 1-j].plot(np.arange(np.min(d_had[1, s]), np.max(d_had[1, s]),0.01),np.arange(np.min(d_had[1, s]), np.max(d_had[1, s]),0.01),color='k')

    rows = ['Winter Hemisphere', 'Summer Hemisphere']
    for ax, row in zip(axs[:, 0], rows):
        ax.annotate(row, xy=(0, 0.5), xytext=(-40, 0),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='right', va='center', rotation=90, fontsize=12)
    for ax, lab in zip(axs[0], ['a', 'b']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=18, fontweight='bold')
    for ax, lab in zip(axs[1], ['c', 'd']):
        ax.annotate(lab, xy=(0, 1), xytext=(6, 3),
                    xycoords='axes fraction', textcoords='offset points',
                    ha='center', va='baseline', fontsize=18, fontweight='bold')


def scatter_sub_vt(d_had, d_jet):
    fig, axs = plt.subplots(1, 2, figsize=(11, 11/2),
                            facecolor='w', edgecolor='k', linewidth=2)
    for s in range(2):
        a_s = stats.linregress(np.ndarray.flatten(d_had[s,:,17:]), np.ndarray.flatten(d_jet[s, :,17:]))
        a_w = stats.linregress(np.ndarray.flatten(d_had[s,:,:17]), np.ndarray.flatten(d_jet[s, :,:17]))
        axs[s].scatter(d_had[s,:,:17],d_jet[s, :,:17],color='b')
        axs[s].scatter(d_had[s,:,17:],d_jet[s, :,17:],color='r')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_w[0]+a_w[1],color='b')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_s[0]+a_s[1],color='r')
        axs[s].set_ylim([np.min(d_jet[s]), np.max(d_jet[s])])
        axs[s].set_xlim([np.min(d_had[s]), np.max(d_had[s])])
        axs[s].text(np.max(d_had[s])-.35*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+.05*abs(np.ptp(d_jet[s])),'Summer: $\mathregular{r^2}=$ %.2f' % a_s[2]**2)
        axs[s].text(np.max(d_had[s])-0.32*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+0.01*abs(np.ptp(d_jet[s])),'Winter: $\mathregular{r^2}=$ %.2f' % a_w[2]**2)
        axs[s].set_xlabel('ITCZ shift (deg)')
        if s == 0:
            axs[s].set_title('Winter Hemisphere')
            #axs[s].set_ylabel('Change in heat flux across Hadley cell extent (ms$\mathregular{^{-1}}$K$\mathregular{^{-1}}$)')
            axs[s].set_ylabel("Change in Hadley cell strength (kg s$^{-1}$)")
        if s == 1:
            axs[s].set_title('Summer Hemisphere')


def scatter_sub3(d_had, d_jet):
    fig, axs = plt.subplots(1, 2, figsize=(11, 11/2),
                            facecolor='w', edgecolor='k', linewidth=2)
    for s in range(2):
        a_s = stats.linregress(np.ndarray.flatten(d_had[s,:,21:]), np.ndarray.flatten(d_jet[s, :,21:]))
        a_h = stats.linregress(np.ndarray.flatten(d_had[s,:,11:21]), np.ndarray.flatten(d_jet[s, :,11:21]))
        a_w = stats.linregress(np.ndarray.flatten(d_had[s,:,:11]), np.ndarray.flatten(d_jet[s, :,:11]))
        axs[s].scatter(d_had[s,:,:11],d_jet[s, :,:11],color='b')
        axs[s].scatter(d_had[s,:,11:21],d_jet[s, :,11:21],color='purple')
        axs[s].scatter(d_had[s,:,21:],d_jet[s, :,21:],color='r')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_w[0]+a_w[1],color='b')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_h[0]+a_h[1],color='purple')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_s[0]+a_s[1],color='r')
        axs[s].set_ylim([np.min(d_jet[s]), np.max(d_jet[s])])
        axs[s].set_xlim([np.min(d_had[s]), np.max(d_had[s])])
        axs[s].text(np.max(d_had[s])-.35*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+.09*abs(np.ptp(d_jet[s])),'Summer: $\mathregular{r^2}=$ %.2f' % a_s[2]**2)
        axs[s].text(np.max(d_had[s])-0.32*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+0.05*abs(np.ptp(d_jet[s])),'Hadley: $\mathregular{r^2}=$ %.2f' % a_h[2]**2)
        axs[s].text(np.max(d_had[s])-0.32*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+0.01*abs(np.ptp(d_jet[s])),'Winter: $\mathregular{r^2}=$ %.2f' % a_w[2]**2)
        axs[s].set_xlabel('ITCZ shift (deg)')
        if s == 0:
            axs[s].set_title('Winter Hemisphere')
            #axs[s].set_ylabel('Change in heat flux across Hadley cell extent (ms$\mathregular{^{-1}}$K$\mathregular{^{-1}}$)')
            axs[s].set_ylabel("Change in Hadley cell strength (kg s$^{-1}$)")
        if s == 1:
            axs[s].set_title('Summer Hemisphere')


def scatter_sub4(d_had, d_jet):
    fig, axs = plt.subplots(1, 2, figsize=(11, 11/2),
                            facecolor='w', edgecolor='k', linewidth=2)
    for s in range(2):
        a_s = stats.linregress(np.ndarray.flatten(d_had[s,:,25:]), np.ndarray.flatten(d_jet[s, :,25:]))
        a_h = stats.linregress(np.ndarray.flatten(d_had[s,:,17:21]), np.ndarray.flatten(d_jet[s, :,17:21]))
        a_h1 = stats.linregress(np.ndarray.flatten(d_had[s,:,21:25]), np.ndarray.flatten(d_jet[s, :,21:25]))
        a_w = stats.linregress(np.ndarray.flatten(d_had[s,:,:17]), np.ndarray.flatten(d_jet[s, :,:17]))
        axs[s].scatter(d_had[s,:,:17],d_jet[s, :,:17],color='b')
        axs[s].scatter(d_had[s,:,17:21],d_jet[s, :,17:21],color='purple')
        axs[s].scatter(d_had[s,:,21:25],d_jet[s, :,21:25],color='green')
        axs[s].scatter(d_had[s,:,25:],d_jet[s, :,25:],color='r')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_w[0]+a_w[1],color='b')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_h[0]+a_h[1],color='purple')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_h[0]+a_h1[1],color='green')
        axs[s].plot(np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01),np.arange(np.min(d_had[s]), np.max(d_had[s]),0.01)*a_s[0]+a_s[1],color='r')
        axs[s].set_ylim([np.min(d_jet[s]), np.max(d_jet[s])])
        axs[s].set_xlim([np.min(d_had[s]), np.max(d_had[s])])
        axs[s].text(np.max(d_had[s])-.35*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+.13*abs(np.ptp(d_jet[s])),'Summer: $\mathregular{r^2}=$ %.2f' % a_s[2]**2)
        axs[s].text(np.max(d_had[s])-0.32*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+0.09*abs(np.ptp(d_jet[s])),'Hadley: $\mathregular{r^2}=$ %.2f' % a_h[2]**2)
        axs[s].text(np.max(d_had[s])-0.32*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+0.05*abs(np.ptp(d_jet[s])),'Hadley1: $\mathregular{r^2}=$ %.2f' % a_h[2]**2)
        
        axs[s].text(np.max(d_had[s])-0.32*abs(np.ptp(d_had[s])),
                        np.min(d_jet[s])+0.01*abs(np.ptp(d_jet[s])),'Winter: $\mathregular{r^2}=$ %.2f' % a_w[2]**2)
        axs[s].set_xlabel('ITCZ shift (deg)')
        if s == 0:
            axs[s].set_title('Winter Hemisphere')
            #axs[s].set_ylabel('Change in heat flux across Hadley cell extent (ms$\mathregular{^{-1}}$K$\mathregular{^{-1}}$)')
            axs[s].set_ylabel("Change in Hadley cell strength (kg s$^{-1}$)")
        if s == 1:
            axs[s].set_title('Summer Hemisphere')


def zonalplot(plotdata, strm, sigma, lat, maxx=1, title=''):
    """
    Plots input grid of quantity at sigma and lat coords

    Parameters
    ----------
    plotdata: array
        data being plotted
    sigma: array
        sigma levels of data
    lat: array
        latitudes of data
    title: str
        optional title
    """
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    plt.figure(facecolor='white')
    #  plt.suptitle("Tropical + polar heating", size=36, fontweight='bold', y=.95)
    meshlat, meshsigma = np.meshgrid(lat, sigma)
    mycmap2 = plt.cm.YlOrRd(np.arange(256))
    mycmap1 = plt.cm.Blues_r(np.arange(256))
    my_cmap = np.concatenate((mycmap1, mycmap2), axis=0)
    my_cmap[230:282, :] = 1
    newcmap = mpl.colors.LinearSegmentedColormap.from_list("newjet", my_cmap)
    ctrs = np.linspace(-maxx, maxx, 17)
    plot = plt.contourf(meshlat, meshsigma, plotdata, ctrs, extend='both',
                        cmap=newcmap, vmin=-maxx, vmax=maxx)
    c = plt.colorbar(plot, orientation='horizontal',
                     shrink=1, spacing='proportional', aspect=50,
                     ticks=np.linspace(-maxx, maxx, 9))
    #plt.contour(meshlat, meshsigma, strm, np.arange(-2, 2, .1), colors='k')
    c.set_label('Potential temperature difference', fontsize='20')
    c.ax.set_xticklabels(np.linspace(-maxx, maxx, 9), fontsize='20')
    plt.gca().invert_yaxis()
    plt.yscale('linear')
    ax = plt.gca()
    ax.set_xlim([-90, 90])
    plt.xlabel('Latitude (deg)', fontsize=20)
    plt.ylabel('Pressure (hPa)', fontsize=20)
    plt.xticks(np.arange(-90, 105, 15), fontsize='20')
    plt.yticks(fontsize='20')
    plt.title(title, y=1.05, fontsize=30)
    plt.show()


def jetindicesdaily(u850, lat):
    """
    Outputs maximum uwnd speed and position of maximum of jet for
    iGCM output

    Parameters
    ----------
    u850: array
        data to determine jet indices for
    lat: array
        latitude values for data

    Returns
    -------
    jetspeed: float
        max jet speed
    jetlatitude: float
        position of maximum
    """
    # from sympy import Symbol
    # from sympy import solve
    from scipy import optimize
    speeds = np.zeros((2, np.ma.size(u850, axis=0)))
    lats = np.zeros((2, np.ma.size(u850, axis=0)))
    for t in range(np.ma.size(u850, axis=0)):
        wd, ww, sd, sw = jetdist(u850[t, :], lat)
        wm, sm = np.argmax(ww), np.argmax(sw)
        wp = np.polyfit(wd[wm-1:wm+2], ww[wm-1:wm+2], 2)
        sp = np.polyfit(sd[sm-1:sm+2], sw[sm-1:sm+2], 2)
        '''
        x = Symbol('x', real=True)
        f = wp[0] * x**2 + wp[1] * x + wp[2]
        fprime = f.diff(x)
        xmax = solve(fprime, x)
        ymax = wp[0] * xmax[0]**2 + wp[1] * xmax[0] + wp[2]
        a = Symbol('a', real=True)
        b = sp[0] * a**2 + sp[1] * a + sp[2]
        bprime = b.diff(a)
        amax = solve(bprime, a)
        bmax = sp[0] * amax[0]**2 + sp[1] * amax[0] + sp[2]
        '''
        def f(x): return wp[0] * x**2 + wp[1] * x + wp[2]
        xmax = optimize.fmin(lambda x: -f(x), wm, disp=False)
        ymax = f(xmax)

        def g(a): return sp[0] * a**2 + sp[1] * a + sp[2]
        amax = optimize.fmin(lambda a: -g(a), sm, disp=False)
        bmax = g(amax)
        speeds[0, t] = ymax
        lats[0, t] = xmax[0]
        speeds[1, t] = bmax
        lats[1, t] = amax[0]
    return speeds, lats


def jetdist(u850, lat):
    """
    Determines the winter and summer jet distributions

    Parameters
    ----------
    u850: array
        wind data at u850
    lat: array
        latitudes of data

    Returns
    -------
    wdist: array
        winter lat distribution
    wweights: array
        speeds at each lat
    sdist: array
        summer lat distribution
    sweights: array
        speeds at each lat
    """
    h = len(lat) // 2
    wwmx = np.max(u850[:h])
    windex = np.argmax(u850[:h])
    wdmax = lat[windex]
    wdist = [wdmax]
    wweights = [wwmx]

    windex -= 1
    while u850[windex+1] >= u850[windex] and u850[windex+1] >= 0 and windex >= 0:
        if u850[windex] < 0:
            lat_0 = lat[windex] - u850[windex] * ((lat[windex+1]-lat[windex]) /
                                                  (u850[windex+1] -
                                                   u850[windex]))
            wdist.insert(0, lat_0)
            wweights.insert(0, 0)
            windex -= 1
        else:
            wdist.insert(0, lat[windex])
            wweights.insert(0, u850[windex])
            windex -= 1

    windex = np.argmax(u850[:h])
    windex += 1
    while u850[windex-1] >= u850[windex] and u850[windex-1] >= 0:
        if u850[windex] < 0:
            lat_0 = lat[windex-1] - u850[windex-1] * ((lat[windex] -
                                                       lat[windex-1]) /
                                                      (u850[windex] -
                                                       u850[windex-1]))
            wdist.append(lat_0)
            wweights.append(0)
            windex += 1
        else:
            wdist.append(lat[windex])
            wweights.append(u850[windex])
            windex += 1

    swmx = np.max(u850[h:])
    sindex = np.argmax(u850[h:]) + h
    sdmax = lat[sindex]
    sdist = [sdmax]
    sweights = [swmx]

    sindex -= 1
    while u850[sindex+1] >= u850[sindex] and u850[sindex+1] >= 0:
        if u850[sindex] < 0:
            lat_0 = lat[sindex] - u850[sindex] * ((lat[sindex+1]-lat[sindex]) /
                                                  (u850[sindex+1] -
                                                   u850[sindex]))
            sdist.insert(0, lat_0)
            sweights.insert(0, 0)
            sindex -= 1
        else:
            sdist.insert(0, lat[sindex])
            sweights.insert(0, u850[sindex])
            sindex -= 1

    sindex = np.argmax(u850[h:]) + h
    sindex += 1
    while sindex < 2*h and u850[sindex-1] >= u850[sindex] and u850[sindex-1] >= 0:
        if u850[sindex] < 0:
            lat_0 = lat[sindex-1] - u850[sindex-1] * ((lat[sindex] -
                                                      lat[sindex-1]) /
                                                      (u850[sindex] -
                                                       u850[sindex-1]))
            sdist.append(lat_0)
            sweights.append(0)
            sindex += 1
        else:
            sdist.append(lat[sindex])
            sweights.append(u850[sindex])
            sindex += 1

    return wdist, wweights, sdist, sweights


