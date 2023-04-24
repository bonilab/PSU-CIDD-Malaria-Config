import yaml
import pandas as pd
import numpy as np
from math import log
import copy


stream = open('input_A1_OZFQ.yml', 'r')
data = yaml.load(stream, Loader=yaml.FullLoader)
stream.close()

data['starting_date'] = '2008/1/1'
data['ending_date'] = '2027/1/1'
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

treatment_coverages = {
    0.55: {
            0.06413: 3,
            0.0585: 2,
            0.0538: 1,
            0.0508: 0.5,
            0.0475: 0.1,
        },
    0.75: {
        0.0735: 3,
        0.067: 2,
        0.061: 1,
        0.057: 0.5,
        0.0537: 0.1,
    },
}


# with importatiation
imp = "_imp"
itc = ""

cost_of_resistance_factors = [1,  2,  5, 10]

mda_coverages = [0.6, 0.7, 0.8, 0.9]

#### Main (A) Drug AL
#### B Drug: AL / (OZ+FQ / KAF-LUM)
mda_therapy_ids = {
        15: [5, 'OZFQ'],
        }

mda_stategy= ['AAB', 'ABA', 'ABB']

oz_efficacies = {
       40: {
               'ec50' : 1.37,
           },
       80: {
               'ec50' : 1.02,
           },
        }
params = []

for oz_eff, ec50_dict in oz_efficacies.items():
    for mda_round in number_MDA_rounds:
        for tc,beta_pfprs in treatment_coverages.items():
            for beta,pfpr in beta_pfprs.items():
                for cr in cost_of_resistance_factors:
                    for mda_coverage in mda_coverages:
                        for mda_therapy_id,mda_therapy in mda_therapy_ids.items():
                            for strategy in mda_stategy:

                                new_data = copy.deepcopy(data)
                                new_data['location_db']['p_treatment_for_less_than_5_by_location'] = np.full(number_of_locations, tc+0.05).tolist()
                                new_data['location_db']['p_treatment_for_more_than_5_by_location'] = np.full(number_of_locations, tc-0.05).tolist()
                                new_data['location_db']['beta_by_location'] = np.full(number_of_locations, beta).tolist()

                                for index,event in enumerate(data['events']):
                                    if event['name'] == 'single_round_MDA':
                                        new_data['events'][index]['info'] = data['events'][index]['info'][0:mda_round]
                                        for i in range(mda_round):
                                            new_data['events'][index]['info'][i]['fraction_population_targeted'] = np.full(number_of_locations, mda_coverage).tolist()

                                for index,event in enumerate(data['events']):
                                    if event['name'] == 'change_treatment_coverage':
                                        new_data['events'][index]['info']= []

                                for locus in new_data["genotype_info"]["loci"]:
                                    for allele in locus["alleles"]:
                                        allele["daily_cost_of_resistance"] = allele["daily_cost_of_resistance"] * cr

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

                                new_data['drug_db'][8]['EC50']['....1'] = ec50_dict['ec50']
                                new_data['drug_db'][9]['EC50']['...1.'] = ec50_dict['ec50']

                                output_filename = 'A3_OZFQ_EF40_80/%d.yml'%(len(params))
                                output_stream = open(output_filename, 'w')
                                yaml.dump(new_data, output_stream)
                                output_stream.close()
                                params.append((mda_round,tc,beta,pfpr,cr,mda_coverage,strategy,mda_therapy[1], oz_eff))
params_df = pd.DataFrame(params, columns=['mda_round','tc','beta','pfpr0','cr_factor','mda_coverage','strategy','b_drug', 'oz_eff'])
params_df.to_csv('A3_OZFQ_EF40_80_params.csv', index=True)
