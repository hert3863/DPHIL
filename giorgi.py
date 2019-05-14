# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 10:52:29 2017

@author: bakerh
"""


'''
outline = dict()
outline["AUS"] = (( 110, -45), ( 155, -45), ( 155, -11), ( 110, -11))
outline["AMZ"] = (( -82, -20), ( -34, -20), ( -34,  12), ( -82,  12))
outline["SSA"] = (( -76, -56), ( -40, -56), ( -40, -20), ( -76, -20))
outline["CAM"] = ((-116,  10), ( -83,  10), ( -83,  25), ( -85,  25), ( -85,  30), (-116,  30))
outline["WNA"] = ((-130,  30), (-103,  30), (-103,  60), (-130,  60))
outline["CNA"] = ((-103,  30), ( -85,  30), ( -85,  50), (-103,  50))
outline["ENA"] = (( -85,  25), ( -60,  25), ( -60,  50), ( -85,  50))
outline["ALA"] = ((-170,  60), (-103,  60), (-103,  72), (-170,  72))
outline["GRL"] = ((-103,  50), ( -10,  50), ( -10,  85), (-103,  85))
outline["MED"] = (( -10,  30), (  40,  30), (  40,  48), ( -10,  48))
outline["NEU"] = (( -10,  48), (  40,  48), (  40,  75), ( -10,  75))
outline["WAF"] = (( -20, -12), (  22, -12), (  22,  18), ( -20,  18))
outline["EAF"] = ((  22, -12), (  52, -12), (  52,  18), (  22,  18))
outline["SAF"] = (( -10, -35), (  52, -35), (  52, -12), ( -10, -12))
outline["SAH"] = (( -20,  18), (  65,  18), (  65,  30), ( -20,  30))
outline["SEA"] = ((  95, -11), ( 155, -11), ( 155,  20), ( 100,  20), ( 100,   5), (  95,   5))
outline["EAS"] = (( 100,  20), ( 145,  20), ( 145,  50), ( 100,  50))
outline["SAS"] = ((  65,   5), ( 100,   5), ( 100,  30), (  65,  30))
outline["CAS"] = ((  40,  30), (  75,  30), (  75,  50), (  40,  50))
outline["TIB"] = ((  75,  30), ( 100,  30), ( 100,  50), (  75,  50))
outline["NAS"] = ((  40,  50), ( 180,  50), ( 180,  70), (  40,  70))
'''


outlinecpdn = dict()
outlinecpdn["AUS"] = [81, 109, 59, 84]
outlinecpdn["AMZ"] = [62, 89, 148, 175]
outlinecpdn["SSA"] = [88, 118, 151, 172]
outlinecpdn["CAM"] = [48, 65, 130, 149]
outlinecpdn["WNA"] = [24, 49, 123, 138]
outlinecpdn["CNA"] = [32, 49, 137, 148]
outlinecpdn["ENA"] = [32, 53, 147, 161]
outlinecpdn["ALA"] = [14, 25, 102, 138]
outlinecpdn["GRL"] = [4, 33, 137, 188]
outlinecpdn["MED"] = [0, 34, 49, 187, 23]  # 0 indicates goes over meridian
outlinecpdn["NEU"] = [0, 12, 35, 187, 23]  # 0 indicates goes over meridian
outlinecpdn["WAF"] = [0, 58, 83, 181, 13]  # 0 indicates goes over meridian
outlinecpdn["EAF"] = [58, 83, 12, 29]
outlinecpdn["SAF"] = [0, 82, 101, 187, 36]  # 0 indicates goes over meridian
outlinecpdn["SAH"] = [0, 48, 59, 181, 36]  # 0 indicates goes over meridian
outlinecpdn["SEA"] = [56, 82, 51, 84]
outlinecpdn["EAS"] = [32, 57, 53, 78]
outlinecpdn["SAS"] = [48, 69, 35, 54]
outlinecpdn["CAS"] = [32, 49, 21, 41]
outlinecpdn["TIB"] = [32, 49, 40, 54]
outlinecpdn["NAS"] = [16, 33, 21, 97]

'''
coefs={}

for i, item in enumerate(outlinecpdn):
    coefs[item]=mlrcluster(tmeans, wbgt95,msk=outlinecpdn[item])
coef_names={}
coef_names['AUS']='Australia'
coef_names['AMZ']='Amazon Basin'
coef_names['SSA']='Southern South America'
coef_names['CAM']='Central America'
coef_names['WNA']='Western North America'
coef_names['CNA']='Central North America'
coef_names['ENA']='Eastern North America'
coef_names['ALA']='Alaska'
coef_names['GRL']='Greenland'
coef_names['MED']='Mediterranean Basin'
coef_names['NEU']='Northern Europe'
coef_names['WAF']='Western Africa'
coef_names['EAF']='Eastern Africa'
coef_names['SAF']='Southern Africa'
coef_names['SAH']='Sahara'
coef_names['SEA']='South East Asia'
coef_names['EAS']='East Asia'
coef_names['SAS']='South Asia'
coef_names['CAS']='Central Asia'
coef_names['NAS']='North Asia'
coef_names['TIB']='Tibet'
'''
