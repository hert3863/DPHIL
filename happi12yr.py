# -*- coding: utf-8 -*-
"""
Created on Tue Feb  7 16:38:39 2017

@author: bakerh
"""

import numpy as np


def genperts(include_zero_pert=False):
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


def genxmlbias1(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlbias1', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>ALLclim_ancil_134months_OSTIA_sst_1985-12-01_1997-01-30</file_sst>\n\
                <file_sice>ALLclim_ancil_134months_OSTIA_ice_1985-12-01_1997-01-30</file_sice>\n\
                <file_so2dms>so2dms_N96_1985_1997</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_N96_1985_1997</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>11</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>1985</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlbias2(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlbias2', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>ALLclim_ancil_134months_OSTIA_sst_1994-12-01_2006-01-30</file_sst>\n\
                <file_sice>ALLclim_ancil_134months_OSTIA_ice_1994-12-01_2006-01-30</file_sice>\n\
                <file_so2dms>so2dms_N96_1994_2006</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_N96_1994_2006</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>11</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>1994</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlpresent(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlpresent', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>ALLclim_ancil_146months_OSTIA_sst_2004-12-01_2017-01-30</file_sst>\n\
                <file_sice>ALLclim_ancil_146months_OSTIA_ice_2004-12-01_2017-01-30</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2004_2017v2</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2004_2017v2</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2004</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlhappi(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlhappi15', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>HAPPI15C_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30</file_sst>\n\
                <file_sice>HAPPI15C_ancil_146months_NEIL_ice_MMM_2088-12-01_2101-01-30</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2088_2101v2</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2088_2101v2</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>CRED_ghg_RCP26_v3_2105</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2088</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlhappilower(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlhappi15l', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>HAPPI15C_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30</file_sst>\n\
                <file_sice>HAPPI15C_ancil_146months_NEIL_ice_MMM_2088-12-01_2101-01-30</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2088_2101v2</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2088_2101v2</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_lower</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2088</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlhappiupper(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlhappi15u', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>HAPPI15C_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30</file_sst>\n\
                <file_sice>HAPPI15C_ancil_146months_NEIL_ice_MMM_2088-12-01_2101-01-30</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2088_2101v2</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2088_2101v2</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_upper</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2088</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlnatural(e_size):
    '''
    Use nat ssts here
    '''
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlnat', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>NATclim_ancil_146months_OSTIA_sst_MMM_2004-12-01_2017-01-30</file_sst>\n\
                <file_sice>NATclim_ancil_14months_OSTIA_ice_0000P</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2004</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlnatsst(e_size):
    '''
    Use nat SSTs here
    '''
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlnatsst', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>NATclim_ancil_146months_OSTIA_sst_MMM_2004-12-01_2017-01-30</file_sst>\n\
                <file_sice>NATclim_ancil_14months_OSTIA_ice_0000P</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2004_2017v2</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2004_2017v2</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2004</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlnatghg(e_size):
    '''
    Use present day SSTs and SICE here
    '''
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlnatghg', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>ALLclim_ancil_146months_OSTIA_sst_2004-12-01_2017-01-30</file_sst>\n\
                <file_sice>ALLclim_ancil_146months_OSTIA_ice_2004-12-01_2017-01-30</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2004</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo


def genxmlhappi2(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlhappi2', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>HAPPI20C_ancil_146months_OSTIA_sst_MMM_2088-12-01_2101-01-30</file_sst>\n\
                <file_sice>HAPPI20C_ancil_146months_NEIL_ice_MMM_2088-12-01_2101-01-30</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2088_2101v2</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2088_2101v2</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_HAPPI2p0_v2</file_ghg>\n\
                <run_years>12</run_years>\n\
                <run_months>1</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2088</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo
