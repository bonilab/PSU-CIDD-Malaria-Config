# -*- coding: utf-8 -*-
"""
Created on Tue Aug  7 14:23:33 2018

@author: NguyenTran
"""

import yaml
import numpy as np
import pandas as pd
import copy

# import inflect;
# p = inflect.engine();


def kFormatter(num):
    return str(num) if num <= 999 else str(round(num/1000)) + 'k'


stream = open('../input_FLAL_v2.yml', 'r')
data = yaml.load(stream, Loader=yaml.FullLoader)
stream.close()

data['starting_date'] = '2008/1/1'
data['ending_date'] = '2042/1/1'
data['start_of_comparison_period'] = '2022/1/1'

data['seasonal_info']['enable'] = 'false'

# 1 location
location_info = [[0, 0, 0]]
number_of_locations = len(location_info)
data['location_db']['location_info'] = location_info


# population size
popsize = 26850
data['location_db']['population_size_by_location'] = [popsize]

# 3RMDA
number_MDA_round = [0, 1, 2, 3, 4]

#
sd_prob_individual_present_at_mda = [0.3, 0.3, 0.3]
data['sd_prob_individual_present_at_mda'] = sd_prob_individual_present_at_mda

# for index,event in enumerate(data['events']):
#    if event['name'] == 'single_round_MDA':
#        data['events'][index]['info'] = data['events'][index]['info'][0:number_MDA_round]


betas = [0.077]

pfpr_beta_map = {0.077: 5}

improved_tc = {True: '_itc',
               False: ''}

importations = {True: '_imp',
                False: ''}

params = []

for mda_round in number_MDA_round:
    for beta in betas:
        for _, itc in improved_tc.items():
            for _, imp in importations.items():
                new_data = copy.deepcopy(data)
                new_data['location_db']['beta_by_location'] = np.full(
                    number_of_locations, beta).tolist()

                for index, event in enumerate(data['events']):
                    if event['name'] == 'single_round_MDA':
                        new_data['events'][index]['info'] = data['events'][index]['info'][0:mda_round]
                pfpr = pfpr_beta_map[beta]

                if itc == '':
                    for index, event in enumerate(data['events']):
                        if event['name'] == 'change_treatment_coverage':
                            new_data['events'][index]['info'] = []
                if imp == '':
                    for index, event in enumerate(data['events']):
                        if event['name'] == 'introduce_parasites_periodically':
                            new_data['events'][index]['info'] = []

                output_filename = 'FigureS2/%d.yml' % (len(params))
                output_stream = open(output_filename, 'w')
                yaml.dump(new_data, output_stream)
                output_stream.close()
                params.append((mda_round, pfpr, itc, imp))

params_df = pd.DataFrame(params, columns=['mda_round', 'pfpr', 'itc', 'imp'])
params_df.to_csv('FigureS2.csv', index=True)
