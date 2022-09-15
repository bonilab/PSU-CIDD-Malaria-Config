# -*- coding: utf-8 -*-
"""
Created on Sept 14 2022

@author: NguyenTran
"""

import yaml;
import numpy as np;
from math import log;
import copy;
import pandas as pd;

#import inflect;
#p = inflect.engine();

def kFormatter(num):
    return str(num) if num <=999 else str(round(num/1000)) +'k';


tcs = {
       0.25: {
           'tc_str': '0p25',
           'pfprs': {
               0.35: 'PFPR25P0',
               0.13: 'PFPR10p0',
               0.049: 'PFPR01p0',
               0.044: 'PFPR00p1'
               }           
           },
       0.5: {
           'tc_str': '0p50',
           'pfprs': {
               0.37: 'PFPR25P0',
               0.15: 'PFPR10p0',
               0.057: 'PFPR01p0',
               0.051: 'PFPR00p1'
               }           
           },
       0.75: {
           'tc_str': '0p75',
           'pfprs': {
               0.42: 'PFPR25P0',
               0.16: 'PFPR10p0',
               0.066: 'PFPR01p0',
               0.0585: 'PFPR00p1'
               }           
           },
       };

firstline_drugs = {
        0: 'DHA-PPQ',
        1: 'ASAQ',
        2: 'AL',
        };

rapt_drugs = {
        7: 'AL',
        8: 'ASAQ',
        9: 'DHA-PPQ'
        };

rapt_periods = [12, 6 , 3];

compliences =[ 0.2, 0.4, 0.6, 0.8, 1.0];

age_start = [5, 18];

params=[];

for tc, tc_details in tcs.items():
    tc_str = tc_details['tc_str'];
    pfprs = tc_details['pfprs'];
    for beta,pfpr in pfprs.items():
        for firstline_drug, firstline_drug_str in firstline_drugs.items():
            for rapt_drug, rapt_drug_str in rapt_drugs.items():
                for rapt_period in rapt_periods:
                    for complience in compliences:
                        for age in age_start:
                            params.append([beta,tc,firstline_drug,rapt_drug, rapt_period, complience, age]);
                            
                            
df_params = pd.DataFrame (params, columns = ['beta','tc', 'firstline_drug', 'rapt_drug', 'rapt_period', 'complience', 'age' ]);
df_params.to_csv("all_params.csv",index=False);

#%%
stream = open('base_input.yml', 'r');
data = yaml.full_load(stream);
stream.close();

run_index_from = 0
run_index_to = 1

#%%
df_params = pd.read_csv("all_params.csv");
for r_index in range(run_index_from, run_index_to):
    row = df_params.iloc[r_index]
    new_data = copy.deepcopy(data)
    new_data['location_db']['beta_by_location'] = np.full(1, row.beta.item()).tolist()
   
    new_data['location_db']['p_treatment_for_less_than_5_by_location'] = np.full(1, row.tc.item()).tolist()
    new_data['location_db']['p_treatment_for_more_than_5_by_location'] = np.full(1, row.tc.item()).tolist()
   
    new_data['strategy_db'][4]['strategy_ids'] = [round(row.firstline_drug.item()), 3]
    new_data['rapt_config']['therapy_id'] = round(row.rapt_drug.item())
    new_data['rapt_config']['period'] = round(row.rapt_period.item())
    new_data['rapt_config']['compliance'] = row.complience.item()
    new_data['rapt_config']['age_start'] = round(row.age.item())
   
    output_filename = 'A1/%d.yml'%(r_index);
    output_stream = open(output_filename, 'w');
    yaml.dump(new_data, output_stream); 
    output_stream.close();
                     
#%%
# A1
filtered_df = df_params[(df_params.tc==0.5)&(df_params.firstline_drug==2)&(df_params.age==5)&(df_params.rapt_drug==9)]
filtered_df.to_csv("A1_params.csv",index=True);
for index, row in filtered_df.iterrows():
    
    new_data = copy.deepcopy(data)
    new_data['location_db']['beta_by_location'] = np.full(1, row.beta.item()).tolist()
   
    new_data['location_db']['p_treatment_for_less_than_5_by_location'] = np.full(1, row.tc.item()).tolist()
    new_data['location_db']['p_treatment_for_more_than_5_by_location'] = np.full(1, row.tc.item()).tolist()
   
    new_data['strategy_db'][4]['strategy_ids'] = [round(row.firstline_drug.item()), 3]
    new_data['rapt_config']['therapy_id'] = round(row.rapt_drug.item())
    new_data['rapt_config']['period'] = round(row.rapt_period.item())
    new_data['rapt_config']['compliance'] = row.complience.item()
    new_data['rapt_config']['age_start'] = round(row.age.item())
   
    output_filename = 'A1/%d.yml'%(index);
    output_stream = open(output_filename, 'w');
    yaml.dump(new_data, output_stream); 
    output_stream.close();
    print('%d.yml'%(index))