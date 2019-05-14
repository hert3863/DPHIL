# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 13:44:23 2017

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


def genxmlpresent(e_size):
    from ANC import ANC
    anc = ANC()
    anc.Start('a001')
    pert_list = genperts()
    wuinfo = {}
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2005_2017</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2005_2017</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2005</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2005_2017</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2005_2017</file_ozone>\n\
                <file_pert>' + pert_list[i+e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2009</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2005_2017</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2005_2017</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2013</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>CRED_ghg_RCP26_v3_2105</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2089</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>CRED_ghg_RCP26_v3_2105</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2093</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>CRED_ghg_RCP26_v3_2105</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2097</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_lower</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2089</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_lower</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2093</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_lower</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2097</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_upper</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2089</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_upper</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2093</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_RCP26_upper</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2097</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>OSICE_natural_0000P</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2005</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>OSICE_natural_0000P</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i+e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2009</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>OSICE_natural_0000P</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2013</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>OSICE_natural_0000P</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2005_2017</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2005_2017</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2005</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>OSICE_natural_0000P</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2005_2017</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2005_2017</file_ozone>\n\
                <file_pert>' + pert_list[i+e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2009</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>OSICE_natural_0000P</file_sice>\n\
                <file_so2dms>so2dms_rcp45_N96_2005_2017</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp45_N96_2005_2017</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>ghg_defaults</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2013</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2005</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i+e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2009</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_prei_N96_1855_0000P</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_preind_N96_1879_0000Pv5</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_1985_2020</file_solar>\n\
                <file_volcanic>volc_1985_2020</file_volcanic>\n\
                <file_ghg>file_ghg_1850</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2013</model_start_year>\n\
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
    runfile = open('/home/bakerh/Documents/DPhil/CPDN/xmlwus', 'w')
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_HAPPI2p0_v2</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2089</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_HAPPI2p0_v2</file_ghg>\n\
                <run_years>5</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2093</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    for i in range(1, e_size+1):
        runfile.write('\
        <experiment>\n\
            <parameters>\n\
                <file_atmos>xhjlya.start.0000.360.new</file_atmos>\n\
                <file_sst>!</file_sst>\n\
                <file_sice>!</file_sice>\n\
                <file_so2dms>so2dms_rcp26_N96_2089_2101</file_so2dms>\n\
                <file_sulphox>oxi.addfa</file_sulphox>\n\
                <file_ozone>ozone_rcp26_N96_2089_2101</file_ozone>\n\
                <file_pert>' + pert_list[i+2*e_size] + '</file_pert>\n\
                <file_solar>solar_marius_2065_2105v2</file_solar>\n\
                <file_volcanic>volc_marius_2065_2105v2</file_volcanic>\n\
                <file_ghg>ghg_HAPPI2p0_v2</file_ghg>\n\
                <run_years>4</run_years>\n\
                <run_months>0</run_months>\n\
                <exptid>' + str(anc.Get()) + '</exptid>\n\
                <model_start_month>12</model_start_month>\n\
                <model_start_year>2097</model_start_year>\n\
            </parameters>\n\
        </experiment>\n')
        wuinfo[anc.Get()] = ('Member ' + str(i), pert_list[i])
        anc.Next()
    runfile.close()
    return wuinfo
