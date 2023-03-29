import yaml
import pandas as pd
import numpy as np
from math import log
import copy


stream = open('input_A1.yml', 'r')
data = yaml.load(stream, Loader=yaml.FullLoader)
stream.close()

data['starting_date'] = '2008/1/1'
data['ending_date'] = '2042/1/1'
data['start_of_comparison_period']= '2022/1/1'

data['seasonal_info']['enable'] = 'false'


# 1 location
location_info =  [[0, 0, 0]]
number_of_locations = len(location_info)
data['location_db']['location_info']= location_info

# 201368

# population size
popsize = 26850
data['location_db']['population_size_by_location'] = [popsize]


# MDA rounds
number_MDA_rounds = [0,1,2,3,4]
#
sd_prob_individual_present_at_mda = [0.3, 0.3, 0.3]
data['sd_prob_individual_present_at_mda'] = sd_prob_individual_present_at_mda

pfprs = {
        0.06413: 'PFPR3p0',
        0.0585: 'PFPR2p0',
        0.0538: 'PFPR1p0',
        0.0508: 'PFPR0p5',
        0.0475: 'PFPR0p1',
        }

# with importatiation
imp = "_imp"
itc = ""

cost_of_resistance_factors = [1,  2,  5, 10]

mda_coverages = [0.6, 0.7, 0.8, 0.9]

#### Main (A) Drug AL
#### B Drug: AL / (OZ+FQ / KAF-LUM)
mda_therapy_ids = {
        9: [3, 'DHAPPQ'],
        12: [4, 'KAFLUM'],
        }

mda_stategy= ['AAB', 'ABA', 'ABB']

params = []

for mda_round in number_MDA_rounds:
    for beta,pfpr_str in pfprs.items():
        for cr in cost_of_resistance_factors:
            for mda_coverage in mda_coverages:
                for mda_therapy_id,mda_therapy in mda_therapy_ids.items():
                    for strategy in mda_stategy:

                        new_data = copy.deepcopy(data)
                        new_data['location_db']['beta_by_location'] = np.full(number_of_locations, beta).tolist()

                        for index,event in enumerate(data['events']):
                            if event['name'] == 'single_round_MDA':
                                new_data['events'][index]['info'] = data['events'][index]['info'][0:mda_round]

                        for index,event in enumerate(data['events']):
                            if event['name'] == 'change_treatment_coverage':
                                new_data['events'][index]['info']= []

                        for locus in new_data["genotype_info"]["loci"]:
                            for allele in locus["alleles"]:
                                allele["daily_cost_of_resistance"] = allele["daily_cost_of_resistance"] * cr

                        new_data["mean_prob_individual_present_at_mda"] = [mda_coverage + 0.05, mda_coverage - 0.05, mda_coverage + 0.05]
                        new_data["sd_prob_individual_present_at_mda"] = [0.2,0.2,0.2]

                        for index,event in enumerate(data['events']):
                            if event['name'] == 'modify_nested_mft_strategy':
                                if strategy == 'AAB':
                                    new_data['mda_therapy_id'] = 7
                                    new_data['events'][index]['info'][0]['strategy_id'] = mda_therapy[0]
                                else:
                                    if strategy == 'ABA':
                                        new_data['mda_therapy_id'] = mda_therapy_id
                                        new_data['events'][index]['info'][0]['strategy_id'] = 1
                                    else:
                                        new_data['mda_therapy_id'] = mda_therapy_id
                                        new_data['events'][index]['info'][0]['strategy_id'] = mda_therapy[0]
                        output_filename = 'A2_20230328/%d.yml'%(len(params))
                        output_stream = open(output_filename, 'w')
                        yaml.dump(new_data, output_stream)
                        output_stream.close()
                        params.append((mda_round,beta,cr,mda_coverage,strategy,mda_therapy[1]))
params_df = pd.DataFrame(params, columns=['mda_round','beta','cr_factor','mda_coverage','strategy','b_drug'])
params_df.to_csv('A2_20230328_params.csv', index=True)
