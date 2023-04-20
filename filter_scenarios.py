"""Filter scenarios from the AR6 database as basis of cross-sector pathway."""
import os

import numpy
import pandas

# directory containing scenario data downloaded from IIASA
_PROJ_DIR = "C:/Users/ginger.kowal/Documents/Scenario review"


def replicate_sr15_filter():
	"""Replicate Andres's filter of SR15 database."""
	sr15_key_path = os.path.join(
		_PROJ_DIR, 'IPCC_SR15/sr15_metadata_indicators_r2.0.xlsx')
	scen_key = pandas.read_excel(sr15_key_path, sheet_name='meta')

	sr15_scen_path = os.path.join(
		_PROJ_DIR, 'IPCC_SR15', 'iamc15_scenario_data_world_r2.0.xlsx')
	scen_em = pandas.read_excel(sr15_scen_path, sheet_name='data')

	# interquartile range: CO2 from low/no overshoot scenarios
	# identify SR15 low/no overshoot scenarios
	sr15_categories = ['Lower 1.5C low overshoot', 'Below 1.5C']
	filtered_scen = scen_key.loc[
		(scen_key['category'].isin(sr15_categories)) &
		(scen_key['Kyoto-GHG|2010 (SAR)'] == 'in range')]

	# summarize CO2 in low/no overshoots, excluding emissions from FLAG,
	# landfill waste, and fluorinated gases
	emissions_subs = scen_em.loc[
		scen_em['Variable'].startswith('Emissions')]
	print('stop here')

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
