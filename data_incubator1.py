#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 19 10:31:51 2018

@author: bakerh
"""
import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.formula.api import ols


# load data assuming NA and Other/NA are counted as missing values
data = pd.read_excel('/home/bakerh/Downloads/ccrb_datatransparencyinitiative_20170207.xlsx', 'Complaints_Allegations', na_values=['NA', 'Other/NA'])
# remove spaces from column names
data.columns = [c.replace(' ', '_') for c in data.columns]
# drop all rows with incomplete data
data_c = data.dropna()

# return UniqueComplaintIds
ids = data_c.UniqueComplaintId.unique()
# count number of UniqueComplaintIds
a1 = len(ids)

# remove duplicated ids
data_cu = data_c[~data_c.UniqueComplaintId.duplicated(keep='first')]
# count number in each borough
count = data_cu.Borough_of_Occurrence.value_counts()
# calculate prop in highest borough
a2 = count.max()/count.sum()

# 2016 pop data
pop_data = {}
pop_data['Manhattan'] = 1643734
pop_data['Brooklyn'] = 2629150
pop_data['Queens'] = 2333054
pop_data['Bronx'] = 1455720
pop_data['Staten Island'] = 476015
pop_data['Outside NYC'] = np.nan
# pop_data['Total'] = 8537673
pdpop_data = pd.Series(pop_data, name='Population')
# select only 2016
data_2016 = data_cu.loc[data_cu['Received_Year'] == 2016]
# count number in each borough
count_2016 = data_2016.Borough_of_Occurrence.value_counts()
# compute per 100k
cper100 = 100000*count_2016/pdpop_data
a3 = cper100.max()

# compute number of years per complaint
ny = data_cu.Close_Year-data_cu.Received_Year
# compute mean number of years
a4 = ny.mean()

# pick out all instances of frisk assuming they are all associated
# with a stop instance too
data_sf = data_c[data_c['Allegation_Description'].str.contains("Frisk")]
# remove non unique entries
data_sfu = data_sf[~data_sf.UniqueComplaintId.duplicated(keep='first')]
# find year of peak
frisk_counts = data_sfu.groupby(['Received_Year']).size()
yrmax = frisk_counts.argmax()
# perform linear regression
regs = stats.linregress(np.arange(yrmax, 2017),
                        frisk_counts.loc[yrmax:2016].values)
a5 = regs[0]*2018+regs[1]

# compute contigency table observed values
tt = len(data_cu.loc[(data_cu['Is_Full_Investigation'] == True) &
                     (data_cu['Complaint_Has_Video_Evidence'] == True)])
tf = len(data_cu.loc[(data_cu['Is_Full_Investigation'] == True) &
                     (data_cu['Complaint_Has_Video_Evidence'] == False)])
ft = len(data_cu.loc[(data_cu['Is_Full_Investigation'] == False) &
                     (data_cu['Complaint_Has_Video_Evidence'] == True)])
ff = len(data_cu.loc[(data_cu['Is_Full_Investigation'] == False) &
                     (data_cu['Complaint_Has_Video_Evidence'] == False)])
# compute contigency table expected values
tt_e = (tt+tf)*(tt+ft)/(tt+tf+ft+ff)
tf_e = (tt+tf)*(tf+ff)/(tt+tf+ft+ff)
ft_e = (ft+ff)*(tt+ft)/(tt+tf+ft+ff)
ff_e = (ft+ff)*(tf+ff)/(tt+tf+ft+ff)
# compute chisquared statistic
chisq = (tt-tt_e)**2/tt_e + (tf-tf_e)**2/tf_e + (ft-ft_e)**2/ft_e + (ff-ff_e)**2/ff_e
a6 = chisq

# number of allegations per id
numb_alleg = data_c.groupby(['UniqueComplaintId']).size().values
# create indicator dataframe
data_c['F'] = np.where(data_c['Allegation_FADO_Type'].str.contains('Force'),
                       1, 0)
data_c['A'] = np.where(data_c['Allegation_FADO_Type'].str.contains('Abuse of Authority'),
                       1, 0)
data_c['D'] = np.where(data_c['Allegation_FADO_Type'].str.contains('Discourtesy'),
                       1, 0)
data_c['O'] = np.where(data_c['Allegation_FADO_Type'].str.contains('Offensive Language'),
                       1, 0)
# create indicator indices
f = np.greater(data_c.groupby(['UniqueComplaintId'])['F'].sum(), 0).astype(int).values
a = np.greater(data_c.groupby(['UniqueComplaintId'])['A'].sum(), 0).astype(int).values
d = np.greater(data_c.groupby(['UniqueComplaintId'])['D'].sum(), 0).astype(int).values
o = np.greater(data_c.groupby(['UniqueComplaintId'])['O'].sum(), 0).astype(int).values
# perform regression forcing through the origin
regdata = pd.DataFrame({'f': f, 'a': a,'d': d,'o': o, 'n_a': numb_alleg})
model = ols("n_a ~ f + a + d + o-1", regdata).fit()
# print(model.summary())
a7 = np.max(model._results.params[1:])

# number of precincts
pr_data = {}
pr_data['Manhattan'] = 22
pr_data['Brooklyn'] = 23
pr_data['Queens'] = 16
pr_data['Bronx'] = 12
pr_data['Staten Island'] = 4
pr_data['Outside NYC'] = np.nan
# pop_data['Total'] = 8537673
pdpr_data = pd.Series(pr_data, name='Precincts')
# compute officers in each borough
n_officers = 36000
n_officers_b = cper100*n_officers/cper100.sum()
n_ofc_p = n_officers_b/pdpr_data
# compute ratio
a8 = n_ofc_p.max()/n_ofc_p.min()

























