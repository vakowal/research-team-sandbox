"""Filter scenarios from the AR6 database as basis of cross-sector pathway."""
import os

import numpy
import pandas

# directory containing scenario data downloaded from IIASA
_PROJ_DIR = "C:/Users/ginger.kowal/Documents/Scenario review"

# GWP100 conversion values from AR5 report
_N2O_GWP100_AR5 = 265
_CH4_GWP100_AR5 = 28


def cross_sector_sr15():
	"""Calculate the cross-sector pathway from SR15 database."""
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

	# filter to low/no overshoot scenarios only
	lno_em = scen_em.loc[scen_em['scen_id'].isin(scen_ids)]

	# filter to get CO2 emissions only
	co2_em = lno_em.loc[
		lno_em['Variable'].str.startswith('Emissions|CO2')]

	# exclude emissions from FLAG, landfill waste, and fluorinated gases
	co2_var = 'Emissions|CO2|Energy and Industrial Processes'
	co2_em = co2_em.loc[co2_em['Variable'] == co2_var]

	# calculate gross CO2 emissions from energy and industrial processes
	# TODO

	# calculate interquartile range of yearly values
	year_col = [col for col in co2_em if col.startswith('2')] + ['Variable']
	sr15_25perc = co2_em[year_col].groupby('Variable').quantile(q=0.25)
	sr15_25perc['scenario_col'] = 'IPCC SR15 (25th percentile)'
	sr15_75perc = co2_em[year_col].groupby('Variable').quantile(q=0.75)
	sr15_75perc['scenario_col'] = 'IPCC SR15 (75th percentile)'

	# add N2O from energy: mean of low/no overshoot scenarios
	n2o_var = 'Emissions|N2O|Energy'
	mean_lno_em = lno_em[year_col].groupby('Variable').mean().reset_index()
	energy_N2O = mean_lno_em.loc[mean_lno_em['Variable'] == n2o_var]

	# convert N2O to CO2e using GWP100 from IPCC AR5
	energy_N2O_CO2eq = energy_N2O * _N2O_GWP100_AR5  # TODO CHECK UNITS

	# add CH4 from NZE
	# convert to CO2e using GWP100 from IPCC AR5

	# add gross CO2, N2O, and CH4: this is the cross-sector pathway


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
	cross_sector_sr15()

if __name__ == '__main__':
    main()
