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

# GW100 conversion values from AR6
_CH4FOSS_GWP100_AR6 = 29.8  # GWP100 for fossil methane

# conversion factor from kt to Mt
_KT_to_MT = 0.001

# conversion factor from Gt to Mt
_GT_to_MT = 1000

# medium concern, yearly deployment of BECCS in 2050 (Gt CO2/year)
_MED_BECCS = 3

# high concern, yearly deployment of BECCS in 2050 (Gt CO2/year)
_HI_BECCS = 7

# medium concern, yearly deployment of fossil CCS in 2050 (Gt CO2/year)
_MED_F_CCS = 3.8

# high concern, yearly deployment of fossil CCS in 2050 (Gt CO2/year)
_HI_F_CCS = 8.8

# maximum yearly sequestration via af-/reforestation in 2050 (Gt CO2/year)
_MAX_AFOLU = 3.6


def fill_EIP_emissions(em_df, id_list):
    """Calculate net energy & industrial process emissions.

    Not all scenarios contain the variable "Emissions|CO2|Energy and Industrial
    Processes". For those that don't, calculate it as the sum of energy, and
    industrial process CO2 emissions.

    Args:
    	em_df (Pandas dataframe): dataframe containing emissions data
    	id_list (list): complete list of scenario ids to fill

	Returns:
		a Pandas dataframe containing all the values in `em_df`, plus
			calculated values for Energy and Industrial Process emissions
			for those scenarios that didn't originally report it

    """
    year_col = [col for col in em_df if col.startswith('2')]
    summary_cols = ['scen_id', 'Variable'] + year_col
    co2_var = 'Emissions|CO2|Energy and Industrial Processes'
    co2_em = em_df.loc[em_df['Variable'] == co2_var][summary_cols]
    em_df.reset_index(inplace=True)

    # not all models report values for this variable: identify
    # those that don't
    estimate_mod = set(id_list).difference(
        set(em_df.loc[em_df['Variable'] == co2_var]['scen_id']))
    print("**note: estimating EIP emissions for {} scenarios".format(
        len(estimate_mod)))

    # for those, calculate the variable as the sum of emissions from Energy,
    # and from Industrial processes
    r_idx = len(em_df)
    for sid in estimate_mod:
        sum_rows = em_df.loc[
            (em_df['scen_id'] == sid) &
            (em_df['Variable'].isin(
                ['Emissions|CO2|Energy',
                'Emissions|CO2|Industrial Processes']))]
        sum_vals = sum_rows[year_col].sum(skipna=False)
        em_df.loc[r_idx] = sum_vals
        em_df.loc[r_idx, 'scen_id'] = sid
        em_df.loc[r_idx, 'Variable'] = co2_var
        r_idx = r_idx + 1

    return em_df


def calc_gross_eip_co2(em_df, year_col):
    """Calculate gross CO2 emissions from EIP.

    Gross CO2 emissions from energy and industrial prcoesses are calculated
    as `Emissions|CO2|Energy and Industrial Processes` +
    `Carbon Sequestration|CCS|Biomass`

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values

    Returns:
        dataframe containing gross EIP CO2 emissions
    """
    ccs_var = 'Carbon Sequestration|CCS|Biomass'
    eip_var = 'Emissions|CO2|Energy and Industrial Processes'
    sum_df = em_df.loc[em_df['Variable'].isin([ccs_var, eip_var])]
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
	sr15_25perc['source'] = 'IPCC SR15 (25th percentile)'
	sr15_25perc['perc'] = '25th perc'
	sr15_50perc = co2_em[quant_col].groupby('Variable').quantile(q=0.5)
	sr15_50perc['source'] = 'IPCC SR15 (median)'
	sr15_50perc['perc'] = 'median'
	sr15_75perc = co2_em[quant_col].groupby('Variable').quantile(q=0.75)
	sr15_75perc['source'] = 'IPCC SR15 (75th percentile)'
	sr15_75perc['perc'] = '75th perc'
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


def read_ar6_data():
	"""Read AR6 scenario data from file.

	Args:
		None

	Returns:
		a tuple, (ar6_key, ar6_scen) containing:
		ar6_key (Pandas dataframe): dataframe containing scenario metadata
		ar6_scen (Pandas dataframe): dataframe containing scenario data

	"""
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

	return ar6_key, ar6_scen


def filter_AR6_scenarios():
	"""Filter scenarios from the AR6 database."""
	ar6_key, ar6_scen = read_ar6_data()
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
	c1_gross_eip_co2 = calc_gross_eip_co2(c1_co2_em, year_col)
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

	iisd_gross_eip_co2 = calc_gross_eip_co2(issd_co2_em, year_col)
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


def calc_ch4_updated():
	"""Calculate updated methane for inclusion in cross-sector pathway."""
	ch4_df = pandas.read_csv(os.path.join(_PROJ_DIR, 'methane summary.csv'))
	ch4_vals = ch4_df.loc[
		ch4_df['Source'] == 'NZE WEO 2022 / methane tracker 2023']
	year_col = [col for col in ch4_vals if col.startswith('2')]
	ch4_co2eq = ch4_vals[year_col] * _CH4FOSS_GWP100_AR6
	# save, or add to CO2 and N20


def extract_imps():
	"""Get net and gross CO2 EIP emissions from the AR6 IMPs."""
	ar6_key, ar6_scen = read_ar6_data()
	year_col = [col for col in ar6_scen if col.startswith('2')]
	C1_IMPs = ['SP', 'LD', 'Ren']
	imp_scen = ar6_key.loc[ar6_key['IMP_marker'].isin(C1_IMPs)]['scen_id']
	imp_em = ar6_scen.loc[ar6_scen['scen_id'].isin(imp_scen)]
	imp_filled = fill_EIP_emissions(imp_em, imp_scen)
	imp_gross = calc_gross_eip_co2(imp_filled, year_col)
	imp_gross.to_csv(
		os.path.join(_OUT_DIR, 'ar6_imp.csv'), index=False)


def implement_filter(scen_id_list, emissions_df, filter_flag):
	"""Filter scenarios according to sustainability/feasibility thresholds.

	Filter a set of scenarios according to 2050 deployment of BECCS, fossil
	CCS, and/or af-/reforestation. Deployment in 2050 is calculated as the
	average of deployment in 2040 and 2060. The specific filters applied are
	determined by the value of `filter_flag`:

		1 (least stringent): remove scenarios exceeding thresholds of high
			concern for 2050 deployment of BECCS and CCS, and exceeding
			sustainable thresholds of af-/reforestation according to the most
			inclusive variable, 'Carbon Sequestration|Land Use'
		2 (medium stringency): remove scenarios exceeding thresholds of medium
			concern for 2050 deployment of BECCS and CCS, and exceeding
			sustainable thresholds of af-/reforestation according to a less
			inclusive variable, 'Carbon Sequestration|Land Use|Afforestation'.
			This approximates the filters of IISD (2021)
		3 (most stringent): remove scenarios exceeding thresholds of medium
			concern for 2050 deployment of BECCS and CCS, and exceeding
			sustainable thresholds of af-/reforestation according to the most
			inclusive variable, 'Carbon Sequestration|Land Use'

	Args:
		scen_id_list (list of strings): list of scenario ids that should be
			filtered. For example, this could be the full list of scenarios
			in the AR6 database, or the subset of C1 scenarios
		emissions_df (Pandas dataframe): dataframe containing data for
			emissions and sequestration, used to identify scenarios meeting
			filtering criteria
		filter_flag (int): flag indicating what filtering criteria should be
			applied

	Returns:
		a list of strings that is a subset of `scen_id_list`, giving the
			scenarios that meet the given filter criteria

	"""
	emissions_df['est2050'] = emissions_df[['2040', '2060']].mean(axis=1)
	if filter_flag == 1:
		beccs_lim = _HI_BECCS
		fccs_lim = _HI_F_CCS
		afolu_var = 'Carbon Sequestration|Land Use'

	elif filter_flag == 2:
		beccs_lim = _MED_BECCS
		fccs_lim = _MED_F_CCS
		afolu_var = 'Carbon Sequestration|Land Use|Afforestation'

	elif filter_flag == 3:
		beccs_lim = _MED_BECCS
		fccs_lim = _MED_F_CCS
		afolu_var = 'Carbon Sequestration|Land Use'

	else:
		raise ValueError("Filter flag must be in {1, 2, 3}")

	rem1 = set(
		emissions_df.loc[
			(emissions_df['Variable'] == 'Carbon Sequestration|CCS|Biomass') &
			(emissions_df['est2050'] > (beccs_lim * _GT_to_MT))]['scen_id'])
	rem2 = set(
		emissions_df.loc[
			(emissions_df['Variable'] == 'Carbon Sequestration|CCS|Fossil') &
			(emissions_df['est2050'] > (fccs_lim * _GT_to_MT))]['scen_id'])
	rem3 = set(
		emissions_df.loc[
			(emissions_df['Variable'] == afolu_var) &
			(emissions_df['est2050'] > (_MAX_AFOLU * _GT_to_MT))]['scen_id'])
	rem_ids = rem1.union(rem2).union(rem3)

	filtered_ids = set(scen_id_list).difference(rem_ids)
	return filtered_ids


def compare_ar6_filters():
	"""Compare different ways of filtering the AR6 scenarios."""
	ar6_key, ar6_scen = read_ar6_data()
	year_col = [col for col in ar6_scen if col.startswith('2')]

	c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
	c1_filter1 = implement_filter(c1_scen, ar6_scen, filter_flag=1)
	c1_filter2 = implement_filter(c1_scen, ar6_scen, filter_flag=2)
	c1_filter3 = implement_filter(c1_scen, ar6_scen, filter_flag=3)

	# fill EIP emissions for all C1 scenarios in the database
	ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)

	# calculate gross emissions for all C1 scenarios in the database
	ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)

	# calculate median gross emissions for scenarios in filtered sets
	summary_cols = ['2020', '2030', '2050', 'Variable']
	med_filter0 = ar6_gross_em.loc[
		ar6_gross_em['scen_id'].isin(c1_scen)][summary_cols].groupby(
			'Variable').quantile(q=0.5)
	med_filter0['filter_flag'] = 0

	med_filter1 = ar6_gross_em.loc[
		ar6_gross_em['scen_id'].isin(c1_filter1)][summary_cols].groupby(
			'Variable').quantile(q=0.5)
	med_filter1['filter_flag'] = 1

	med_filter2 = ar6_gross_em.loc[
		ar6_gross_em['scen_id'].isin(c1_filter2)][summary_cols].groupby(
			'Variable').quantile(q=0.5)
	med_filter2['filter_flag'] = 2

	med_filter3 = ar6_gross_em.loc[
		ar6_gross_em['scen_id'].isin(c1_filter3)][summary_cols].groupby(
			'Variable').quantile(q=0.5)
	med_filter3['filter_flag'] = 3
	f_df = pandas.concat([med_filter0, med_filter1, med_filter2, med_filter3])
	f_df['percch_2030'] = (f_df['2030'] - f_df['2020']) / f_df['2020']
	f_df['percch_2050'] = (f_df['2050'] - f_df['2020']) / f_df['2020']
	f_df.to_csv(
		os.path.join(_OUT_DIR, "ar6_filtered_sets.csv"), index=False)


def main():
	# cross_sector_sr15()
	# filter_AR6_scenarios()
	# calc_ch4_updated()
	# extract_imps()
	compare_ar6_filters()


if __name__ == '__main__':
    main()
