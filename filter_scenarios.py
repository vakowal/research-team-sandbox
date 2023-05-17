"""Filter scenarios from the AR6 database as basis of cross-sector pathway."""
import os

import numpy
import pandas

# directory containing scenario data downloaded from IIASA
_PROJ_DIR = "C:/Users/ginger.kowal/Documents/Scenario review"

# path to mirrored files on google drive
_GDRIVE = "G:/.shortcut-targets-by-id/1rSoiKOBotDMn7VymKwxdpRQDr7ixLAhv"

# output directory
_OUT_DIR = os.path.join(
	_GDRIVE,
	"Research Team/04-Current projects/1.5C Scenarios Review 2023/Intermediate analysis products")

# GWP100 conversion values from AR5 report
_N2O_GWP100_AR5 = 265
_CH4_GWP100_AR5 = 28

# conversion factor from kt to Mt
_KT_to_MT = 0.001


def fill_EIP_emissions(em_df, id_list):
    """Calculate net energy & industrial process emissions."""
    year_col = [col for col in em_df if col.startswith('2')]
    summary_cols = ['scen_id', 'Variable'] + year_col
    co2_var = 'Emissions|CO2|Energy and Industrial Processes'
    co2_em = em_df.loc[em_df['Variable'] == co2_var][summary_cols]
    co2_em.reset_index(inplace=True)

    # not all models report values for this variable: identify
    # those that don't
    estimate_mod = set(id_list).difference(
        set(em_df.loc[em_df['Variable'] == co2_var]['scen_id']))
    print("**note: estimating EIP emissions for {} scenarios".format(
        len(estimate_mod)))

    # for those, calculate the variable as the sum of emissions from Energy,
    # and from Industrial processes
    r_idx = len(co2_em)
    for sid in estimate_mod:
        sum_rows = em_df.loc[
            (em_df['scen_id'] == sid) &
            (em_df['Variable'].isin(
                ['Emissions|CO2|Energy',
                'Emissions|CO2|Industrial Processes']))]
        sum_vals = sum_rows[year_col].sum(skipna=False)
        co2_em.loc[r_idx] = sum_vals
        co2_em.loc[r_idx, 'scen_id'] = sid
        co2_em.loc[r_idx, 'Variable'] = co2_var
        r_idx = r_idx + 1

    return co2_em


def calc_gross_eip_co2(net_co2_df, em_df, year_col):
    """Calculate gross CO2 emissions from EIP.

    Gross CO2 emissions from energy and industrial prcoesses are calculated
    as `Emissions|CO2|Energy and Industrial Processes` +
    `Carbon Sequestration|CCS|Biomass`

    Args:
        net_co2_df (pandas dataframe): dataframe containing
            net Energy and Industrial Process emissions
        em_df (pandas dataframe): dataframe containing all
            other variables

    Returns:
        dataframe containing gross EIP CO2 emissions
    """
    ccs_var = 'Carbon Sequestration|CCS|Biomass'
    ccs_df = em_df.loc[em_df['Variable'] == ccs_var]
    sum_df = pandas.concat([ccs_df, net_co2_df])
    sum_cols = year_col + ['scen_id']
    gross_eip_co2_df = sum_df[sum_cols].groupby('scen_id').sum()
    gross_eip_co2_df.reset_index(inplace=True)
    gross_eip_co2_df['Variable'] = 'Emissions|CO2|Energy and Industrial Processes|Gross'
    return gross_eip_co2_df


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

	# calculate gross CO2 emissions from energy and industrial processes,
	# excluding emissions from FLAG, landfill waste, and fluorinated gases
	year_col = [col for col in lno_em if col.startswith('2')]
	sum_col = year_col + ['scen_id']
	sum_var = [
		'Emissions|CO2|Energy and Industrial Processes',
		'Carbon Sequestration|CCS|Biomass']
	co2_em = lno_em.loc[lno_em['Variable'].isin(sum_var)]
	gross_co2 = co2_em[sum_col].groupby('scen_id').sum()
	gross_co2['Variable'] = 'Emissions|CO2|Energy and Industrial Processes|Gross'
	gross_co2.reset_index(inplace=True)
	co2_em = pandas.concat([co2_em, gross_co2])

	# calculate interquartile range of yearly values
	quant_col = year_col + ['Variable']
	sr15_25perc = co2_em[quant_col].groupby('Variable').quantile(q=0.25)
	sr15_25perc['scenario_col'] = 'IPCC SR15 (25th percentile)'
	sr15_50perc = co2_em[quant_col].groupby('Variable').quantile(q=0.5)
	sr15_50perc['scenario_col'] = 'IPCC SR15 (median)'
	sr15_75perc = co2_em[quant_col].groupby('Variable').quantile(q=0.75)
	sr15_75perc['scenario_col'] = 'IPCC SR15 (75th percentile)'
	sr15_df = pandas.concat(
		[sr15_25perc, sr15_50perc, sr15_75perc]).reset_index()
	sr15_df.to_csv(
		os.path.join(_OUT_DIR, "sr15_filtered_summary.csv"), index=False)

	# add N2O from energy: mean of low/no overshoot scenarios
	# note that this is probably incorrect, Andres likely summed N2O from more
	# than just the Energy sector.
	n2o_var = 'Emissions|N2O|Energy'
	mean_lno_em = lno_em[quant_col].groupby('Variable').mean().reset_index()
	energy_N2O = mean_lno_em.loc[mean_lno_em['Variable'] == n2o_var]

	# convert N2O to CO2e using GWP100 from IPCC AR5
	energy_N2O_CO2eq = energy_N2O[year_col] * _KT_to_MT * _N2O_GWP100_AR5
	energy_N2O_CO2eq['scenario_col'] = 'IPCC SR15 (mean)'
	energy_N2O_CO2eq['Variable'] = 'Emissions|N2O|Energy_CO2eq'

	# add CH4 from NZE: directly from Andres's notes in PtN-Z supplementary
	# data document. This is taken from fig. 3.5 in NZE or equivalent, and
	# already converted to CO2eq using the GWP100 value from AR5
	ch4_df = pandas.read_csv(
		os.path.join(_PROJ_DIR, 'IEA/IEA_NZE_fossil_methane_Andres.csv'))
	ch4_df['scenario_col'] = 'IEA_NZE'
	ch4_df['Variable'] = 'Emissions|CH4|Fossil_CO2eq'

	# add gross CO2, N2O, and CH4: this is the cross-sector pathway
	cross_sector_df = pandas.concat([sr15_df, energy_N2O_CO2eq, ch4_df])
	# cross_sector_df.to_csv(
	#  	"C:/Users/ginger.kowal/Desktop/cross_sector_sr15.csv")


def filter_AR6_scenarios():
	"""Filter scenarios from the AR6 database."""
	key_path = os.path.join(
		_PROJ_DIR, 'IPCC_AR6',
		'AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx')
	ar6_key = pandas.read_excel(
		key_path, sheet_name='meta_Ch3vetted_withclimate')
	ar6_key['scen_id'] = ar6_key['Model'] + ' ' + ar6_key['Scenario']

	scen_path = os.path.join(
		_PROJ_DIR, 'IPCC_AR6/AR6_Scenarios_Database_World_v1.1.csv')
	ar6_scen = pandas.read_csv(scen_path)
	ar6_scen['scen_id'] = ar6_scen['Model'] + ' ' + ar6_scen['Scenario']
	year_col = [col for col in ar6_scen if col.startswith('2')]

	# select all C1 scenarios: 1.5C with low or no overshoot
	c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
	c1_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_scen)]

	# calculate net CO2 emissions from energy and industrial processes
	c1_co2_em = fill_EIP_emissions(c1_em, c1_scen)

	# Filter C1 scenarios according to feasibility, following IISD 2022:
	# Exclude scenarios exceeding the “medium” feasibility thresholds for the
	# scaling up of fossil fuel CCS, bioenergy with CCS (BECCS), and
	# afforestation or reforestation.
	iisd_scen_path = os.path.join(
	    _PROJ_DIR, 'IISD/IISD_2022_filtered_AR6_scenarios.csv')
	iisd_scen = pandas.read_csv(iisd_scen_path)
	iisd_scen_ids = iisd_scen['Model Scenario']

	iisd_em = ar6_scen.loc[ar6_scen['scen_id'].isin(iisd_scen_ids)]
	issd_co2_em = fill_EIP_emissions(iisd_em, iisd_scen_ids)

	# compare net CO2 emissions for C1 and IISD filtered scenarios
	summary_cols = year_col + ['Variable']

	c1_25perc = c1_co2_em[summary_cols].groupby('Variable').quantile(q=0.25)
	c1_25perc['source'] = 'AR6 C1'
	c1_25perc['perc'] = '25th perc'
	c1_75perc = c1_co2_em[summary_cols].groupby('Variable').quantile(q=0.75)
	c1_75perc['source'] = 'AR6 C1'
	c1_75perc['perc'] = '75th perc'
	c1_med = c1_co2_em[summary_cols].groupby('Variable').quantile(q=0.5)
	c1_med['source'] = 'AR6 C1'
	c1_med['perc'] = 'median'

	iisd_25perc = issd_co2_em[summary_cols].groupby('Variable').quantile(q=0.25)
	iisd_25perc['source'] = 'AR6 C1+IISD'
	iisd_25perc['perc'] = '25th perc'
	iisd_75perc = issd_co2_em[summary_cols].groupby('Variable').quantile(q=0.75)
	iisd_75perc['source'] = 'AR6 C1+IISD'
	iisd_75perc['perc'] = '75th perc'
	iisd_med = issd_co2_em[summary_cols].groupby('Variable').quantile(q=0.5)
	iisd_med['source'] = 'AR6 C1+IISD'
	iisd_med['perc'] = 'median'

	# Calculate gross CO2 emissions from energy and industrial processes
	c1_gross_eip_co2 = calc_gross_eip_co2(c1_co2_em, c1_em, year_col)
	c1_gr_25perc = c1_gross_eip_co2[
		summary_cols].groupby('Variable').quantile(q=0.25)
	c1_gr_25perc['source'] = 'AR6 C1'
	c1_gr_25perc['perc'] = '25th perc'
	c1_gr_75perc = c1_gross_eip_co2[
		summary_cols].groupby('Variable').quantile(q=0.75)
	c1_gr_75perc['source'] = 'AR6 C1'
	c1_gr_75perc['perc'] = '75th perc'
	c1_gr_med = c1_gross_eip_co2[
		summary_cols].groupby('Variable').quantile(q=0.5)
	c1_gr_med['source'] = 'AR6 C1'
	c1_gr_med['perc'] = 'median'

	iisd_gross_eip_co2 = calc_gross_eip_co2(issd_co2_em, iisd_em, year_col)
	iisd_gr_25perc = iisd_gross_eip_co2[
		summary_cols].groupby('Variable').quantile(q=0.25)
	iisd_gr_25perc['source'] = 'AR6 C1+IISD'
	iisd_gr_25perc['perc'] = '25th perc'
	iisd_gr_75perc = iisd_gross_eip_co2[
		summary_cols].groupby('Variable').quantile(q=0.75)
	iisd_gr_75perc['source'] = 'AR6 C1+IISD'
	iisd_gr_75perc['perc'] = '75th perc'
	iisd_gr_med = iisd_gross_eip_co2[
		summary_cols].groupby('Variable').quantile(q=0.5)
	iisd_gr_med['source'] = 'AR6 C1+IISD'
	iisd_gr_med['perc'] = 'median'

	ar6_df = pandas.concat([
	    c1_25perc, c1_75perc, c1_med, iisd_25perc, iisd_75perc, iisd_med,
	    c1_gr_25perc, c1_gr_75perc, c1_gr_med, iisd_gr_25perc,
	    iisd_gr_75perc, iisd_gr_med])
	ar6_df.reset_index(inplace=True)
	ar6_df.to_csv(
		os.path.join(_OUT_DIR, "ar6_df.csv"), index=False)


def main():
	# cross_sector_sr15()
	filter_AR6_scenarios()


if __name__ == '__main__':
    main()
