#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 10:43:58 2019

@author: bakerh
"""
import numpy as np


def interest1(nheld, yr, ntimes, reinvest='off'):
    bonds = 71787616800
    init = nheld
    intrate = np.zeros(ntimes)
    for j in range(ntimes):
        nheld = init
        total_winnings = 0
        for i in range(12*yr):
            total = 0
            w_np = np.random.rand(1)
            if nheld*2 > w_np * bonds:
                total += 1000000
            if nheld*5 > w_np * bonds:
                total += 100000
            if nheld*11 > w_np * bonds:
                total += 50000
            if nheld*20 > w_np * bonds:
                total += 25000
            if nheld*51 > w_np * bonds:
                total += 10000
            if nheld*103 > w_np * bonds:
                total += 5000
            if nheld*1831 > w_np * bonds:
                total += 1000
            if nheld*5493 > w_np * bonds:
                total += 500
            if nheld*25194 > w_np * bonds:
                total += 100
            if nheld*25194 > w_np * bonds:
                total += 50
            if nheld*3145004 > w_np * bonds:
                total += 25
            total_winnings += total
            if reinvest == 'on':
                nheld = total + nheld
                total_winnings = 0
        intrate[j] = 100*np.exp(np.log((nheld+total_winnings)/init)/yr)-100
    return intrate.mean()


def interest2(nheld, yr, ntimes, reinvest='off'):
    bonds = 71787616800
    init = nheld
    intrate = np.zeros(ntimes)
    for j in range(ntimes):
        nheld = init
        total_winnings = 0
        for i in range(12*yr):
            total = 0
            w_np = np.random.rand(nheld)
            win = np.where(2 > w_np * bonds, 1, 0)
            total += np.sum(win)*1000000
            win = np.where(5 > w_np * bonds, 1, 0)
            total += np.sum(win)*100000
            win = np.where(11 > w_np * bonds, 1, 0)
            total += np.sum(win)*50000
            win = np.where(20 > w_np * bonds, 1, 0)
            total += np.sum(win)*25000
            win = np.where(51 > w_np * bonds, 1, 0)
            total += np.sum(win)*10000
            win = np.where(103 > w_np * bonds, 1, 0)
            total += np.sum(win)*5000
            win = np.where(1831 > w_np * bonds, 1, 0)
            total += np.sum(win)*1000
            win = np.where(5493 > w_np * bonds, 1, 0)
            total += np.sum(win)*500
            win = np.where(25194 > w_np * bonds, 1, 0)
            total += np.sum(win)*100
            win = np.where(25194 > w_np * bonds, 1, 0)
            total += np.sum(win)*50
            win = np.where(3145004 > w_np * bonds, 1, 0)
            total += np.sum(win)*25
            total_winnings += total
            if reinvest == 'on':
                nheld = total + nheld
                total_winnings = 0
        intrate[j] = 100*np.exp(np.log((nheld+total_winnings)/init)/yr)-100
    return intrate.mean()


def gaussian(n):
    import matplotlib.pyplot as plt
    market = np.zeros(n)
    for i in range(len(market)-1):
        market[i+1] = market[i] + np.random.normal()
    plt.plot(market)
    return market