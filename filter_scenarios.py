"""Filter scenarios from the AR6 database as basis of cross-sector pathway."""
import os

import numpy
import pandas

# directory containing scenario data downloaded from IIASA
_PROJ_DIR = "C:/Users/ginger.kowal/Documents/Scenario review"


def replicate_sr15_filter():
	"""Replicate Andres's filter of SR15 database."""
	# identify SR15 low/no overshoot scenarios
	sr15_key_path = os.path.join(
		_PROJ_DIR, 'IPCC_SR15/sr15_metadata_indicators_r2.0.xlsx')
	scen_key = pandas.read_excel(sr15_key_path, sheet_name='meta')

	sr15_categories = ['1.5C low overshoot', 'Below 1.5C']
	filtered_scen = scen_key.loc[
		(scen_key['category'].isin(sr15_categories)) &
		(scen_key['Kyoto-GHG|2010 (SAR)'] == 'in range')]
	scen_ids = filtered_scen['model'] + filtered_scen['scenario']

	# CO2 in low/no overshoot scenarios
	sr15_scen_path = os.path.join(
		_PROJ_DIR, 'IPCC_SR15', 'iamc15_scenario_data_world_r2.0_data.csv')
	scen_em = pandas.read_csv(sr15_scen_path)
	scen_em['scen_id'] = scen_em['Model'] + scen_em['Scenario']

	# filter to get CO2 emissions only
	filtered_em = scen_em.loc[
		scen_em['Variable'].str.startswith('Emissions|CO2')]

	# filter to low/no overshoot scenarios only
	filtered_em = filtered_em.loc[
		filtered_em['scen_id'].isin(scen_ids)]

	# calculate interquartile range
	sr15_25perc = filtered_em.groupby('Variable').quantile(q=0.25)
	sr15_25perc['id'] = 'IPCC SR15 (25th percentile)'
	sr15_75perc = filtered_em.groupby('Variable').quantile(q=0.75)
	sr15_75perc['id'] = 'IPCC SR15 (75th percentile)'

	# exclude emissions from FLAG, landfill waste, and fluorinated gases

	# add N20 from energy: mean of low/no overshoot scenarios
	# convert to CO2e using GWP100 from IPCC AR5

	# add CH4 from NZE
	# convert to CO2e using GWP100 from IPCC AR5


def filter_AR6_scenarios():
	"""Filter scenarios from the AR6 database."""
	key_path = os.path.join(
		_PROJ_DIR, 'IPCC_AR6',
		'AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx')
	ar6_key = pandas.read_excel(
		key_path, sheet_name='meta_Ch3vetted_withclimate')

	scen_path = os.path.join(
		_PROJ_DIR, 'IPCC_AR6/AR6_Scenarios_Database_World_v1.1.csv')
	ar6_scen = pandas.read_csv(scen_path)


def main():
	replicate_sr15_filter()

if __name__ == '__main__':
    main()
