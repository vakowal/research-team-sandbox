"""Filter scenarios from the AR6 database as basis of cross-sector pathway."""
import os

import numpy
import pandas

import aneris
from aneris.tutorial import load_data


# directory containing scenario data downloaded from IIASA
_PROJ_DIR = "C:/Users/Ginger.Kowal.CORPCDP/Documents/Scenario review"

# output directory
_OUT_DIR = os.path.join(_PROJ_DIR, 'outputs')


# GWP100 conversion values from AR5 report
_N2O_GWP100_AR5 = 265
_CH4_GWP100_AR5 = 28

# GW100 conversion values from AR6
_CH4FOSS_GWP100_AR6 = 29.8  # GWP100 for fossil methane
_N2O_GWP100_AR6 = 273
_HFC_GWP100_AR6 = 1530
_PFC_GWP100_AR6 = 7380
_SF6_GWP100_AR6 = 25200

# cumulative emissions from FLAG, 2020-2050 (GtCO2e)
_FLAG_2020_50_CO2E = -99.54

# conversion factor from kt to Mt
_KT_to_MT = 0.001

# conversion factor from Gt to Mt
_GT_to_MT = 1000

# medium concern, yearly energy from bioenergy in any year between 2010-2050
# (EJ/year)
_MED_BIO = 100.

# medium concern, yearly deployment of BECCS in 2050 (Gt CO2/year)
_MED_BECCS = 3

# high concern, yearly deployment of BECCS in 2050 (Gt CO2/year)
_HI_BECCS = 7

# medium concern, yearly deployment of fossil CCS in 2050 (Gt CO2/year)
_MED_F_CCS = 3.8

# high concern, yearly deployment of fossil CCS in 2050 (Gt CO2/year)
_HI_F_CCS = 8.8

# maximum yearly sequestration via af-/reforestation in any year between
# 2010-2050 (Gt CO2/year)
_MAX_AFOLU = 3.6

# maximum cumulative CCS between 2010 and 2050 (Gt)
_MAX_CCS = 214

# gross energy and industrial process CO2 emissions in 2020 (Mt CO2)
# (Forster et al 2023)
_2020_CO2 = 35289.751

# maximum CDR via novel methods in 2020 (Mt CO2) (State of CDR Report)
_2020_CDR = 2.3

# year of historical emissions to be used for harmonization
_HARM_YEAR = 2022


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
        sum_vals = sum_rows[year_col].sum()
        em_df.loc[r_idx] = sum_vals
        em_df.loc[r_idx, 'scen_id'] = sid
        em_df.loc[r_idx, 'Variable'] = co2_var
        r_idx = r_idx + 1

    return em_df


def calc_gross_eip_co2(em_df, year_col, fill_df=False):
    """Calculate gross CO2 emissions from EIP.

    Gross CO2 emissions from energy and industrial prcoesses are calculated
    as `Emissions|CO2|Energy and Industrial Processes` +
    `Carbon Sequestration|CCS|Biomass` +
    `Carbon Sequestration|Direct Air Capture` +
    `Carbon Sequestration|Enhanced Weathering`

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values
        fill_df (bool): fill NaNs in em_df via interpolation?

    Returns:
        dataframe containing gross EIP CO2 emissions

    """
    ccs_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    # check for negative values in CCS. if any values are negative,
    # change the sign
    ccs_df = em_df.loc[em_df['Variable'].isin(ccs_var_list)]
    if (ccs_df[year_col] < 0).any().any():
        for year in year_col:
            ccs_df.loc[ccs_df[year] < 0, year] = -(
                ccs_df.loc[ccs_df[year] < 0, year])
    eip_var = 'Emissions|CO2|Energy and Industrial Processes'
    eip_df = em_df.loc[em_df['Variable'] == eip_var]
    sum_df = pandas.concat([ccs_df, eip_df])
    sum_df.replace(0, numpy.nan, inplace=True)
    if fill_df:
        filled_df = sum_df[year_col].interpolate(axis=1)
    else:
        filled_df = sum_df
    filled_df['scen_id'] = sum_df['scen_id']
    sum_cols = year_col + ['scen_id']

    gross_eip_co2_df = filled_df[sum_cols].groupby('scen_id').sum()
    gross_eip_co2_df.reset_index(inplace=True)
    gross_eip_co2_df['Variable'] = 'Emissions|CO2|Energy and Industrial Processes|Gross'
    return gross_eip_co2_df


def calc_eip_n2o(em_df, year_col):
    """Calculate N2O emissions from energy and industrial processes.

    Energy-related (i.e., non-AFOLU) N2O emissions are calculated as:
    'Emissions|N2O|Energy' + 'Emissions|N2O|Industrial Processes' +
    'Emissions|N2O|Other' + 'Emissions|N2O|Waste'.

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values

    Returns:
        A pandas dataframe containing energy-related N2O emissions for each
            scenario, in kt N2O per year

    """
    n2o_var_list = [
        'Emissions|N2O|Energy', 'Emissions|N2O|Industrial Processes',
        'Emissions|N2O|Other', 'Emissions|N2O|Waste']
    sum_cols = year_col + ['scen_id']
    n2o_df = em_df.loc[em_df['Variable'].isin(n2o_var_list)]
    eip_n2o_df = n2o_df[sum_cols].groupby('scen_id').sum()
    return eip_n2o_df


def calc_eip_ch4(em_df, year_col):
    """Calculate CH4 emissions from energy and industrial processes.

    Energy-related (i.e., non-AFOLU) CH4 emissions are calculated as:
    'Emissions|CH4|Energy' + 'Emissions|CH4|Industrial Processes' +
    'Emissions|CH4|Other' + 'Emissions|N2O|Waste'.

    Args:
        em_df (pandas dataframe): dataframe containing emissions variables
        year_col (list): list of columns giving yearly values

    Returns:
        A pandas dataframe containing energy and industrial process CH4
            emissions for each scenario, in Mt CH4 per year

    """
    ch4_var_list = [
        'Emissions|CH4|Energy', 'Emissions|CH4|Industrial Processes',
        'Emissions|CH4|Other', 'Emissions|CH4|Waste']
    sum_cols = year_col + ['scen_id']
    ch4_df = em_df.loc[em_df['Variable'].isin(ch4_var_list)]
    eip_ch4_df = ch4_df[sum_cols].groupby('scen_id').sum()
    return eip_ch4_df


def harmonize_to_historical(
        emissions_df, hist_path, regions_path, config_path):
    """Use aneris to harmonize emissions from AR6 scenarios to 2022 emissions.

    Use the budget conservation method to derive harmonized CO2 emissions
    trajectories from the modeled trajectories in `emissions_df` to 2022
    emissions, conserving the same cumulative emissions budget. The emissions
    budget is calculated over the time period 2022-2100.

    Args:
        emissions_df (pandas dataframe): dataframe of modeled EIP CO2 emissions
        hist_path (path): path to csv file containing historical data
        regions_path (path): path to csv file containing aneris regions and
            sector data
        config_path (path): path to yaml file containing config info for aneris

    Returns:
        a dictionary containing the following structures:
            dict['harmonized'] harmonized emissions trajectories
            dict['metadata'] aneris metadata
            dict['diagnostics'] aneris diagnostics

    """
    # last year of harmonization period: defines the time period for cumulative
    # emissions that should be maintained during harmonization
    max_year = 2051
    year_col = [str(i) for i in range(2020, max_year)]
    an_mod_cols = ([
        'Model', 'Scenario', 'Region', 'Variable', 'Unit', 'scen_id'] +
        year_col)

    # Interpolate yearly modeled data
    raw_mod = emissions_df[an_mod_cols]
    raw_mod.replace(0, numpy.nan, inplace=True)
    modeled_data = raw_mod[year_col].interpolate(axis=1)
    harm_cols = [str(i) for i in range(_HARM_YEAR, max_year)]
    modeled_data = modeled_data[harm_cols]
    modeled_data['Model'] = 'model'
    modeled_data['Region'] = 'World'
    modeled_data['Unit'] = 'Mt CO2/yr'
    modeled_data['Variable'] = (
        'p|Emissions|CO2|Energy and Industrial Processes|s')
    modeled_data['Scenario'] = raw_mod['scen_id']

    # historical emissions and aneris config files
    hist_data = pandas.read_csv(hist_path)
    aneris_regions = pandas.read_csv(regions_path)
    aneris_rc = aneris.RunControl(rc=config_path)
    aneris_rc['config']['harmonize_year'] = _HARM_YEAR

    overrides = pandas.DataFrame(
        [], columns=['Model', 'Scenarimodo', 'Region', 'Variable', 'Unit'])
    aneris_driver = aneris.HarmonizationDriver(
        aneris_rc, hist_data, modeled_data, overrides, aneris_regions)
    # carbon budget conservation
    aneris_driver.overrides = modeled_data[
        ["Model", "Scenario", "Region", "Variable", "Unit"]].assign(
            Method="budget")

    # harmonize
    failed_list = []
    er_list = []
    for scenario in aneris_driver.scenarios():
        try:
            aneris_driver.harmonize(scenario)
        except ValueError as error:
            failed_list.append(scenario)
            er_list.append(error)

    harmonized, metadata, diagnostics = aneris_driver.harmonized_results()
    ret_dict = {
        'harmonized': harmonized,
        'metadata': metadata,
        'diagnostics': diagnostics,
    }
    return ret_dict


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
    #   "C:/Users/ginger.kowal/Desktop/cross_sector_sr15.csv")


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


def NZE_eip_ch4_co2eq():
    """Calculate updated methane for inclusion in cross-sector pathway.

    Read a csv file containing energy-related methane emissions from several
    sources. Select methane from the NZE WEO 2022 and convert it to Mt CO2eq
    using the conversion factor for fossil methane from AR6 report.

    Returns:
        A pandas dataframe containing energy-related methane emissions for the
            years 2020, 2021, 2025, 2030, 2040, and 2050, in Mt CO2eq per year

    """
    ch4_df = pandas.read_csv(os.path.join(_OUT_DIR, 'methane summary.csv'))
    ch4_vals = ch4_df.loc[
        ch4_df['Source'] == 'NZE WEO 2022 / methane tracker 2023']
    year_col = [col for col in ch4_vals if col.startswith('2')]
    ch4_co2eq = ch4_vals[year_col] * _CH4FOSS_GWP100_AR6
    return ch4_co2eq


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


def calc_afolu_co2(emissions_df):
    id_df = emissions_df.set_index('scen_id')
    test_col = [str(idx) for idx in list(range(2020, 2051))]
    # calculate total AFOLU CO2e
    sum_cols = test_col + ['scen_id']
    co2_df = id_df.loc[id_df['Variable'] == 'Emissions|CO2|AFOLU']
    co2_df.replace(0, numpy.nan, inplace=True)
    afolu_co2 = co2_df[test_col].interpolate(axis=1).sum(axis=1)
    afolu_co2_df = pandas.DataFrame(afolu_co2)
    afolu_co2_df.reset_index(inplace=True)
    return afolu_co2_df


def calc_afolu_co2e(emissions_df, ch4_gwp100, n2o_gwp100):
    """Calculate AFOLU emissions in CO2eq in 2050.

    Args:
        emissions_df (pandas dataframe): dataframe containing emissions
            data
        ch4_gwp100 (float): GWP100 value to use for CH4
        n2o_gwp100 (float): GWP100 value to use for N2O

    Returns:
        a pandas series containing total AFOLU emissions in 2050, in MtCO2e

    """
    test_col = [str(idx) for idx in list(range(2020, 2051))]
    # calculate total AFOLU CO2e
    co2_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CO2|AFOLU']
    ch4_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CH4|AFOLU']
    ch4_co2eq = ch4_df[test_col] * ch4_gwp100
    n2o_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|N2O|AFOLU'][
            test_col].groupby('scen_id').sum()
    n2o_co2eq = n2o_df * n2o_gwp100 * _KT_to_MT
    n2o_co2eq.reset_index(inplace=True)

    co2e_df = pandas.concat(
        [co2_df, ch4_df, n2o_df]).groupby('scen_id').sum()
    afolu_2050_co2e = co2e_df['2050']
    return afolu_2050_co2e


def flag_filter(emissions_df):
    """Filter scenarios for compatibility with SBTi FLAG pathway.

    The SBTi FLAG pathway uses mitigation potentials from Roe et al (2019) and
    a baseline emissions value (from Roe et al) to calculate a 1.5C-compatible
    emissions pathway for FLAG. Use this pathway to identify scenarios where
    emissions from the land sector are smaller than this in terms of cumulative
    CO2e between 2020 and 2050.  CH4 and N2O are converted to CO2eq using
    GWP-100 values from AR5 (comparable to the FLAG pathway).

    Args:
        emissions_df (Pandas dataframe): dataframe containing data for
            emissions and sequestration, used to identify scenarios meeting
            filtering criteria

    Returns:
        a list of strings, giving the scenarios that should be removed
            according to compatibility with the FLAG pathway

    """
    test_col = [str(idx) for idx in list(range(2020, 2051))]
    # calculate total AFOLU CO2e
    sum_cols = test_col + ['scen_id']
    co2_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CO2|AFOLU']
    ch4_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CH4|AFOLU']
    ch4_co2eq = ch4_df[sum_cols] * _CH4_GWP100_AR5
    n2o_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|N2O|AFOLU'][
            sum_cols].groupby('scen_id').sum()
    n2o_co2eq = n2o_df * _N2O_GWP100_AR5 * _KT_to_MT
    n2o_co2eq.reset_index(inplace=True)

    co2e_df = pandas.concat([co2_df, ch4_df, n2o_df]).groupby('scen_id').sum()
    co2e_df.replace(0, numpy.nan, inplace=True)
    co2e_df.reset_index(inplace=True)

    afolu_em = co2e_df[test_col].interpolate(axis=1).sum(axis=1)
    afolu_limit = _FLAG_2020_50_CO2E  * _GT_to_MT
    rem_afolu = set(
        co2e_df.loc[co2e_df[test_col].interpolate(
            axis=1).sum(axis=1) < afolu_limit]['scen_id'])

    return rem_afolu


def sustainability_filters(scen_id_list, emissions_df, filter_flag):
    """Filter scenarios according to first draft sustainability thresholds.

    Filter a set of scenarios according to 2050 deployment of biofuels, CCS,
    and af-/reforestation. Deployment in 2050 is calculated as the average of
    deployment in 2040 and 2060.

    Args:
        scen_id_list (list of strings): list of scenario ids that should be
            filtered. For example, this could be the full list of scenarios
            in the AR6 database, or the subset of C1 scenarios
        emissions_df (Pandas dataframe): dataframe containing data for
            emissions and sequestration, used to identify scenarios meeting
            filtering criteria
        filter_flag (int): flag indicating what combination of filters to apply

    Returns:
        a list of strings that is a subset of `scen_id_list`, giving the
            scenarios that meet the given filter criteria

    """
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    biom_df = emissions_df.loc[
        emissions_df['Variable'] == 'Primary Energy|Biomass']
    rem_biom = set(
        biom_df.loc[(biom_df[test_col] > _MED_BIO).any(axis=1)]['scen_id'])

    # lu_var = 'Carbon Sequestration|Land Use'
    lu_var = 'Carbon Sequestration|Land Use|Afforestation'
    lu_df = emissions_df.loc[emissions_df['Variable'] == lu_var]
    rem_lu = set(
        lu_df.loc[
            (lu_df[test_col] > (_MAX_AFOLU * _GT_to_MT)).any(axis=1)][
            'scen_id'])

    ccs_df = emissions_df.loc[
        emissions_df['Variable'] == 'Carbon Sequestration|CCS']
    rem_ccs = set(
        ccs_df.loc[ccs_df[test_col].interpolate(
            axis=1).sum(axis=1) > (_MAX_CCS * _GT_to_MT)]['scen_id'])

    rem_afolu = flag_filter(emissions_df)

    # filter scenarios with >2.3 Mt CDR in 2020
    ccs_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    cdr_df = emissions_df.loc[emissions_df['Variable'].isin(ccs_var_list)]
    cdr_sum_df = cdr_df.groupby('scen_id').sum()
    cdr_sum_df.reset_index(inplace=True)
    rem_cdr = set(cdr_sum_df.loc[cdr_sum_df['2020'] > _2020_CDR]['scen_id'])

    if filter_flag == 1:
        # remove scenarios according to bioenergy only
        rem_ids = rem_biom

    elif filter_flag == 2:
        # remove scenarios according to land use sequestration only
        rem_ids = rem_lu

    elif filter_flag == 3:
        # remove scenarios according to total CCS only
        rem_ids = rem_ccs

    elif filter_flag == 4:
        # remove scenarios according to bioenergy and CCS
        rem_ids = rem_biom.union(rem_ccs)

    elif filter_flag == 5:
        # remove scenarios according to bioenergy, land use, and total CCS
        rem_ids = rem_biom.union(rem_lu).union(rem_ccs)

    elif filter_flag == 6:
        # remove scenarios according to bioenergy, total CCS, and FLAG pathway
        rem_ids = rem_biom.union(rem_ccs).union(rem_afolu)

    elif filter_flag == 7:
        # remove scenarios according to bioenergy, afforestation, total CCS,
        # CDR in 2020, and FLAG pathway
        rem_ids = rem_biom.union(rem_lu).union(rem_ccs).union(rem_afolu).union(
            rem_cdr)

    filtered_ids = set(scen_id_list).difference(rem_ids)
    return filtered_ids


def iisd_filter_variations(scen_id_list, emissions_df, filter_flag):
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
    c1_filter1 = sustainability_filters(c1_scen, ar6_scen, filter_flag=1)
    c1_filter2 = sustainability_filters(c1_scen, ar6_scen, filter_flag=2)
    c1_filter3 = sustainability_filters(c1_scen, ar6_scen, filter_flag=3)
    c1_filter4 = sustainability_filters(c1_scen, ar6_scen, filter_flag=4)
    c1_filter5 = sustainability_filters(c1_scen, ar6_scen, filter_flag=5)
    c1_filter6 = sustainability_filters(c1_scen, ar6_scen, filter_flag=6)
    c1_filter7 = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # fill EIP emissions for all C1 scenarios in the database
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)

    # calculate gross emissions for all C1 scenarios in the database
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)

    # calculate median gross emissions for scenarios in filtered sets
    summary_cols = ['2020', '2030', '2040', '2050', '2060', 'Variable']
    med_co2_filter0 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter0['filter_flag'] = 0

    med_co2_filter1 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter1)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter1['filter_flag'] = 1

    med_co2_filter2 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter2)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter2['filter_flag'] = 2

    med_co2_filter3 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter3)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter3['filter_flag'] = 3

    med_co2_filter4 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter4)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter4['filter_flag'] = 4

    med_co2_filter5 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter5)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter5['filter_flag'] = 5

    med_co2_filter6 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter6)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter6['filter_flag'] = 6

    med_co2_filter7 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter7)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filter7['filter_flag'] = 7

    co2_df = pandas.concat(
        [med_co2_filter0, med_co2_filter1, med_co2_filter2, med_co2_filter3,
        med_co2_filter4, med_co2_filter5, med_co2_filter6, med_co2_filter7])
    co2_df.reset_index(inplace=True)
    co2_df['percch_2030'] = (co2_df['2030'] - co2_df['2020']) / co2_df['2020']
    co2_df['percch_2050'] = (co2_df['2050'] - co2_df['2020']) / co2_df['2020']
    co2_df['num scenarios'] = [
        len(scen_set) for scen_set in [
            c1_scen, c1_filter1, c1_filter2, c1_filter3, c1_filter4,
            c1_filter5, c1_filter6, c1_filter7]]
    co2_df.to_csv(
        os.path.join(_OUT_DIR, "ar6_filtered_sets_co2.csv"), index=False)

    # add CO2eq from CH4 and N2O to gross CO2
    n2o_co2eq_df = calc_eip_n2o_co2eq(ar6_scen, year_col)

    med_n2o_filter0 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter0['filter_flag'] = 0

    med_n2o_filter1 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter1)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter1['filter_flag'] = 1

    med_n2o_filter2 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter2)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter2['filter_flag'] = 2

    med_n2o_filter3 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter3)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter3['filter_flag'] = 3

    med_n2o_filter4 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter4)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter4['filter_flag'] = 4

    med_n2o_filter5 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter5)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter5['filter_flag'] = 5

    med_n2o_filter6 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter6)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter6['filter_flag'] = 6

    med_n2o_filter7 = n2o_co2eq_df.loc[
        n2o_co2eq_df['scen_id'].isin(c1_filter7)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_n2o_filter7['filter_flag'] = 7

    n2o_df = pandas.concat(
        [med_n2o_filter0, med_n2o_filter1, med_n2o_filter2, med_n2o_filter3,
        med_n2o_filter4, med_n2o_filter5, med_n2o_filter6, med_n2o_filter7])
    n2o_df.reset_index(inplace=True)
    n2o_df.to_csv(
        os.path.join(_OUT_DIR, 'ar6_n2o_filtered_sets.csv'), index=False)

    ch4_co2eq_df = calc_eip_ch4_co2eq(ar6_scen, year_col)

    med_ch4_filter0 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter0['filter_flag'] = 0

    med_ch4_filter1 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter1)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter1['filter_flag'] = 1

    med_ch4_filter2 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter2)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter2['filter_flag'] = 2

    med_ch4_filter3 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter3)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter3['filter_flag'] = 3

    med_ch4_filter4 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter4)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter4['filter_flag'] = 4

    med_ch4_filter5 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter5)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter5['filter_flag'] = 5

    med_ch4_filter6 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter6)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter6['filter_flag'] = 6

    med_ch4_filter7 = ch4_co2eq_df.loc[
        ch4_co2eq_df['scen_id'].isin(c1_filter7)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_ch4_filter7['filter_flag'] = 7

    ch4_df = pandas.concat(
        [med_ch4_filter0, med_ch4_filter1, med_ch4_filter2, med_ch4_filter3,
        med_ch4_filter4, med_ch4_filter5, med_ch4_filter6, med_ch4_filter7])
    ch4_df['Units'] = 'MtCO2e/year'
    ch4_df.reset_index(inplace=True)
    ch4_df.to_csv(
        os.path.join(_OUT_DIR, 'ar6_ch4_filtered_sets.csv'), index=False)

    cs_df = pandas.concat(
        [co2_df, n2o_df, ch4_df]).groupby('filter_flag').sum()
    cs_df.reset_index(inplace=True)
    cs_df['percch_2030'] = (cs_df['2030'] - cs_df['2020']) / cs_df['2020']
    cs_df['percch_2050'] = (cs_df['2050'] - cs_df['2020']) / cs_df['2020']
    cs_df['num scenarios'] = [
        len(scen_set) for scen_set in [
            c1_scen, c1_filter1, c1_filter2, c1_filter3, c1_filter4,
            c1_filter5, c1_filter6, c1_filter7]]
    cs_df.to_csv(
        os.path.join(_OUT_DIR, "ar6_filtered_sets_cs.csv"), index=False)


def export_data_for_fig():
    """Export gross fossil CO2 from all C1, and from filtered scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    summary_years = ['2020', '2030', '2040', '2050', '2060']

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filter7 = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # fill EIP emissions for all C1 scenarios in the database
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)

    # total CO2 for C1 scenarios
    summary_cols = summary_years + ['scen_id']
    c1_co2 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_scen)][summary_cols]
    c1_co2.reset_index(inplace=True)

    # total CO2eq for filtered scenarios
    med_co2_filtered_c1 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filter7)][
            summary_years].quantile(q=0.5)
    med_co2_filtered_c1['scen_id'] = 'Median of filtered scenarios'
    filtered_df = pandas.DataFrame(med_co2_filtered_c1).transpose()

    fig_df = pandas.concat([c1_co2, filtered_df])
    fig_df.to_csv(os.path.join(_OUT_DIR, 'co2_c1_filtered.csv'), index=False)


def afforestation_test():
    """Why does afforestation filter lead to decreased ambition?"""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']

    # fill EIP emissions for all C1 scenarios in the database
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)

    # calculate gross emissions for all C1 scenarios in the database
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)

    c1_em = ar6_filled_em.loc[ar6_filled_em['scen_id'].isin(c1_scen)]
    c1_gross_em = ar6_gross_em.loc[ar6_gross_em['scen_id'].isin(c1_scen)]
    c1_em.set_index(['scen_id'], inplace=True)
    c1_gross_em.set_index(['scen_id'], inplace=True)

    # calculate max sequestration via afforestation between 2010 and 2050
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    lu_df = c1_em.loc[
        c1_em['Variable'] == 'Carbon Sequestration|Land Use']
    max_lu = lu_df[test_col].max(axis=1)

    net_2020 = c1_em.loc[
        c1_em['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes']['2020']
    gross_2020 = c1_gross_em.loc[
        c1_gross_em['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross']['2020']
    net_2030 = c1_em.loc[
        c1_em['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes']['2030']
    gross_2030 = c1_gross_em.loc[
        c1_gross_em['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross']['2030']
    net_2050 = c1_em.loc[
        c1_em['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes']['2050']
    gross_2050 = c1_gross_em.loc[
        c1_gross_em['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross']['2050']
    test_df = pandas.DataFrame({
        'max_lu_seq': max_lu,
        'net_2020_EIP_CO2': net_2020,
        'gross_2020_EIP_CO2': gross_2020,
        'net_2030_EIP_CO2': net_2030,
        'gross_2030_EIP_CO2': gross_2030,
        'net_2050_EIP_CO2': net_2050,
        'gross_2050_EIP_CO2': gross_2050})
    test_df['perc_ch_net_2030'] = (
        (test_df['net_2030_EIP_CO2'] - test_df['net_2020_EIP_CO2']) /
        test_df['net_2020_EIP_CO2'])
    test_df['perc_ch_gross_2030'] = (
        (test_df['gross_2030_EIP_CO2'] - test_df['gross_2020_EIP_CO2']) /
        test_df['gross_2020_EIP_CO2'])
    test_df['perc_ch_net_2050'] = (
        (test_df['net_2050_EIP_CO2'] - test_df['net_2020_EIP_CO2']) /
        test_df['net_2020_EIP_CO2'])
    test_df['perc_ch_gross_2050'] = (
        (test_df['gross_2050_EIP_CO2'] - test_df['gross_2020_EIP_CO2']) /
        test_df['gross_2020_EIP_CO2'])
    test_df.to_csv(
        os.path.join(_OUT_DIR, "max_lu_seq_vs_net_gross_EIP_CO2.csv"))


def summarize_final_energy():
    """Summarize final energy demand in AR6 scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filter7 = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    summary_cols = ['2020', '2030', '2040', '2050', '2060', 'Variable']
    energy_df = ar6_scen.loc[ar6_scen['Variable'] == 'Final Energy']

    energy_c1_25p = energy_df.loc[
        energy_df['scen_id'].isin(c1_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.25)
    energy_c1_25p['filter_flag'] = 0
    energy_c1_25p['quantile'] = 0.25

    energy_c1_50p = energy_df.loc[
        energy_df['scen_id'].isin(c1_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.50)
    energy_c1_50p['filter_flag'] = 0
    energy_c1_50p['quantile'] = 0.50

    energy_c1_75p = energy_df.loc[
        energy_df['scen_id'].isin(c1_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.75)
    energy_c1_75p['filter_flag'] = 0
    energy_c1_75p['quantile'] = 0.75

    energy_f7_25p = energy_df.loc[
        energy_df['scen_id'].isin(c1_filter7)][summary_cols].groupby(
            'Variable').quantile(q=0.25)
    energy_f7_25p['filter_flag'] = 7
    energy_f7_25p['quantile'] = 0.25

    energy_f7_50p = energy_df.loc[
        energy_df['scen_id'].isin(c1_filter7)][summary_cols].groupby(
            'Variable').quantile(q=0.50)
    energy_f7_50p['filter_flag'] = 7
    energy_f7_50p['quantile'] = 0.50

    energy_f7_75p = energy_df.loc[
        energy_df['scen_id'].isin(c1_filter7)][summary_cols].groupby(
            'Variable').quantile(q=0.75)
    energy_f7_75p['filter_flag'] = 7
    energy_f7_75p['quantile'] = 0.75

    energy_sum = pandas.concat(
        [energy_c1_25p, energy_c1_50p, energy_c1_75p, energy_f7_25p,
        energy_f7_50p, energy_f7_75p])
    energy_sum.reset_index(inplace=True)
    energy_sum.to_csv(
        os.path.join(_OUT_DIR, 'final_energy_summary.csv'), index=False)


def summarize_filtered_key_var():
    """Summarize key variables from the median of filtered scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)
    filtered_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_filtered)]
    filtered_em.set_index('scen_id', inplace=True)
    summary_cols = ['2050', 'Variable']

    # Net CO2e AFOLU emissions in 2050
    afolu_co2e = calc_afolu_co2e(filtered_em, _CH4_GWP100_AR5, _N2O_GWP100_AR5)

    # Maximum yearly primary energy from bioenergy, 2010-2050
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    biom_df = filtered_em.loc[
        filtered_em['Variable'] == 'Primary Energy|Biomass']
    max_biom = biom_df[test_col].max(axis=1)

    # Cumulative CCS, 2010-2050
    ccs_df = filtered_em.loc[
        filtered_em['Variable'] == 'Carbon Sequestration|CCS']
    cum_ccs = ccs_df[test_col].interpolate(axis=1).sum(axis=1)

    # Cumulative gross fossil CO2, 2020-2050
    em_col = [str(idx) for idx in list(range(2020, 2051))]
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    filt_gross_em = ar6_gross_em.loc[ar6_gross_em['scen_id'].isin(c1_filtered)]
    filt_gross_em.set_index('scen_id', inplace=True)
    filt_gross_em.replace(0, numpy.nan, inplace=True)
    cum_gross_co2 = filt_gross_em[em_col].interpolate(axis=1).sum(axis=1)

    # Total atmospheric CDR in 2050
    cdr_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    cdr_df = filtered_em.loc[
        filtered_em['Variable'].isin(cdr_var_list)][summary_cols]
    cdr_sum = cdr_df.groupby('scen_id').sum()['2050']

    # Share of primary energy from renewables in 2050
    c1_prien_2050 = filtered_em.loc[
        filtered_em['Variable'] == 'Primary Energy']['2050']
    c1_renen_2050 = filtered_em.loc[
        filtered_em[
            'Variable'] == 'Primary Energy|Renewables (incl. Biomass)']['2050']
    c1_en_df = pandas.DataFrame({
        'Primary Energy': c1_prien_2050,
        'Primary Energy|Renewables': c1_renen_2050})
    ren_share_2050 = (c1_en_df['Primary Energy|Renewables'] /
        c1_en_df['Primary Energy'])

    key_var_vals = pandas.DataFrame({
        'AFOLU CO2e 2050': afolu_co2e,
        'Max yearly bioenergy': max_biom,
        'Cumulative CCS 2010-2050': cum_ccs,
        'CDR 2050': cdr_sum,
        'Ren share 2050': ren_share_2050,
        'Cumulative gross fossil CO2 2020-2050': cum_gross_co2,
        })
    key_var_med = key_var_vals.quantile(q=0.5)
    key_var_df = pandas.DataFrame({
        'median of filtered scenarios': key_var_med})
    key_var_df.to_csv(
        os.path.join(_OUT_DIR, 'AR6_filtered_key_var_summary.csv'))


def summarize_c1_key_var():
    """Summarize key variables in AR6 C1 scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_scen)]
    c1_em.set_index('scen_id', inplace=True)
    summary_cols = ['2050', 'Variable']

    # Net CO2e AFOLU emissions in 2050
    afolu_co2e = calc_afolu_co2e(c1_em, _CH4_GWP100_AR5, _N2O_GWP100_AR5)

    # Maximum yearly primary energy from bioenergy, 2010-2050
    test_col = [str(idx) for idx in list(range(2010, 2051))]
    biom_df = c1_em.loc[c1_em['Variable'] == 'Primary Energy|Biomass']
    max_biom = biom_df[test_col].max(axis=1)

    # Cumulative CCS, 2010-2050
    ccs_df = c1_em.loc[c1_em['Variable'] == 'Carbon Sequestration|CCS']
    cum_ccs = ccs_df[test_col].interpolate(axis=1).sum(axis=1)

    # Cumulative gross fossil CO2, 2020-2050
    em_col = [str(idx) for idx in list(range(2020, 2051))]
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    c1_gross_em = ar6_gross_em.loc[ar6_gross_em['scen_id'].isin(c1_scen)]
    c1_gross_em.set_index('scen_id', inplace=True)
    c1_gross_em.replace(0, numpy.nan, inplace=True)
    cum_gross_co2 = c1_gross_em[em_col].interpolate(axis=1).sum(axis=1)

    # Total atmospheric CDR in 2050
    cdr_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    cdr_df = c1_em.loc[
        c1_em['Variable'].isin(cdr_var_list)][summary_cols]
    cdr_sum = cdr_df.groupby('scen_id').sum()['2050']

    # Share of primary energy from renewables in 2050
    c1_prien_2050 = c1_em.loc[c1_em['Variable'] == 'Primary Energy']['2050']
    c1_renen_2050 = c1_em.loc[
        c1_em['Variable'] == 'Primary Energy|Renewables (incl. Biomass)'][
            '2050']
    c1_en_df = pandas.DataFrame({
        'Primary Energy': c1_prien_2050,
        'Primary Energy|Renewables': c1_renen_2050})
    ren_share_2050 = (c1_en_df['Primary Energy|Renewables'] /
        c1_en_df['Primary Energy'])

    key_var_vals = pandas.DataFrame({
        'AFOLU CO2e 2050': afolu_co2e,
        'Max yearly bioenergy': max_biom,
        'Cumulative CCS 2010-2050': cum_ccs,
        'CDR 2050': cdr_sum,
        'Ren share 2050': ren_share_2050,
        'Cumulative gross fossil CO2 2020-2050': cum_gross_co2,
        })
    key_var_25p = key_var_vals.quantile(q=0.25)
    key_var_75p = key_var_vals.quantile(q=0.75)
    key_var_df = pandas.DataFrame({
        'C1 25 perc': key_var_25p,
        'C1 75 perc': key_var_75p,
        })
    key_var_df.to_csv(
        os.path.join(_OUT_DIR, 'AR6_C1_key_var_summary.csv'))


def afolu_co2e_ngfs():
    """Calculate AFOLU emissions for NGFS Net Zero scenarios."""
    emissions_df = pandas.read_csv(
        os.path.join(_PROJ_DIR, "NGFS/ngfs_afolu_emissions.csv"))

    test_col = [str(idx) for idx in list(range(2020, 2051, 5))]
    # calculate total AFOLU CO2e
    sum_cols = test_col + ['Model']
    co2_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CO2|AFOLU']
    ch4_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|CH4|AFOLU']
    ch4_co2eq = ch4_df[sum_cols] * _CH4_GWP100_AR5
    n2o_df = emissions_df.loc[
        emissions_df['Variable'] == 'Emissions|N2O|AFOLU'][
            sum_cols].groupby('Model').sum()
    n2o_co2eq = n2o_df * _N2O_GWP100_AR5 * _KT_to_MT
    n2o_co2eq.reset_index(inplace=True)

    co2e_df = pandas.concat([co2_df, ch4_df, n2o_df]).groupby('Model').sum()
    co2e_df.replace(0, numpy.nan, inplace=True)

    afolu_em = co2e_df[test_col].interpolate(axis=1).sum(axis=1)
    afolu_co2 = co2_df[test_col].interpolate(axis=1).sum(axis=1)
    print(afolu_em)


def cross_sector_benchmarks_Sept_2023():
    """Calculate cross sector benchmarks from AR6 and key hybrid scenarios.

    This analysis was followed for the draft of the cross-sector pathway
    revision that was shared with the Scientific Advisory Group in September
    2023. It does not include harmonization.

    """
    # Median of filtered scenarios from AR6 database
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # fill EIP emissions for all C1 scenarios in the database
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)

    # calculate median gross emissions for scenarios in filtered sets
    summary_cols = ['2020', '2030', '2040', '2050', 'Variable']
    num_cols = ['2020', '2030', '2040', '2050']
    med_co2_filtered_c1 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filtered)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filtered_c1.reset_index(inplace=True)
    med_co2_filtered_c1['Source'] = 'filtered C1'

    # gross EIP CO2 emissions in IMPs
    focal_imp = ['LD', 'Ren']
    imp_scen = ar6_key.loc[ar6_key['IMP_marker'].isin(focal_imp)]['scen_id']
    imp_em = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(imp_scen)][num_cols]
    imp_em['Source'] = focal_imp

    # gross EIP CO2 emissions in key external scenarios
    key_scenario_df = pandas.read_csv(
        os.path.join(_PROJ_DIR, 'gross_eip_co2_emissions_focal_scen.csv'))
    can_focal_df = pandas.concat([key_scenario_df, imp_em])

    # select subset of focal scenarios
    focal_scen = ['NZE'] + focal_imp
    focal_df = can_focal_df.loc[can_focal_df['Source'].isin(focal_scen)]

    scen_df = pandas.concat([med_co2_filtered_c1, focal_df])
    scen_df['percch_2030'] = (
        scen_df['2030'] - scen_df['2020']) / scen_df['2020']
    scen_df['percch_2040'] = (
        scen_df['2040'] - scen_df['2020']) / scen_df['2020']
    scen_df['percch_2050'] = (
        scen_df['2050'] - scen_df['2020']) / scen_df['2020']
    scen_df['Variable'] = 'Gross fossil CO2'

    # estimate gross EIP CO2 emissions in 2030, 2040, 2050 from average
    # % change among included scenarios
    co2_df = pandas.DataFrame({
        '2020': [_2020_CO2],
        '2030': [_2020_CO2 + (_2020_CO2 * scen_df['percch_2030'].mean())],
        '2040': [_2020_CO2 + (_2020_CO2 * scen_df['percch_2040'].mean())],
        '2050': [_2020_CO2 + (_2020_CO2 * scen_df['percch_2050'].mean())],
        })
    co2_df['Variable'] = 'Gross fossil CO2'

    # summarize single-gas pathways for non-CO2 GHGs
    eip_n2o_df = calc_eip_n2o(ar6_scen, year_col)
    eip_n2o_df.reset_index(inplace=True)
    med_c1_n2o = eip_n2o_df.loc[
        eip_n2o_df['scen_id'].isin(c1_scen)][num_cols].quantile(q=0.5)
    n2o_df = pandas.DataFrame(med_c1_n2o).transpose()
    n2o_df['Variable'] = 'Fossil N2O'
    n2o_df['Source'] = 'Cross sector benchmark'

    eip_ch4_df = calc_eip_ch4(ar6_scen, year_col)
    eip_ch4_df.reset_index(inplace=True)
    med_c1_ch4 = eip_ch4_df.loc[
        eip_ch4_df['scen_id'].isin(c1_scen)][num_cols].quantile(q=0.5)
    ch4_df = pandas.DataFrame(med_c1_ch4).transpose()
    ch4_df['Variable'] = 'Fossil CH4'
    ch4_df['Source'] = 'Cross sector benchmark'

    med_c1_hfc = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|HFC') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)
    med_c1_pfc = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|PFC') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)
    med_c1_sf6 = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|SF6') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)

    # add CO2eq from non-CO2 GHGs
    non_co2_df = pandas.DataFrame(
        {'N2O (Mt CO2e)': med_c1_n2o * _N2O_GWP100_AR6 * _KT_to_MT,
        'CH4 (Mt CO2e)': med_c1_ch4 * _CH4FOSS_GWP100_AR6,
        'HFC (Mt CO2e)': med_c1_hfc * _HFC_GWP100_AR6 * _KT_to_MT,
        'PFC (Mt CO2e)': med_c1_pfc * _PFC_GWP100_AR6 * _KT_to_MT,
        'SF6 (Mt CO2e)': med_c1_sf6 * _SF6_GWP100_AR6 * _KT_to_MT}).transpose()
    non_co2_df['Variable'] = non_co2_df.index
    non_co2_df['Source'] = 'Cross sector benchmark'

    co2e_df = pandas.concat([co2_df, non_co2_df])
    sum_ser = co2e_df[['2020', '2030', '2040', '2050']].sum()
    sum_df = pandas.DataFrame(sum_ser).transpose()
    sum_df['Variable'] = 'Gross fossil CO2e'
    sum_df['Source'] = 'Cross sector benchmark'

    # summary for inclusion in the report
    co2_df['Variable'] = 'Gross fossil CO2'
    co2_df['Source'] = 'Cross sector benchmark'
    summary_df = pandas.concat([
        scen_df[['Variable', 'Source', '2020', '2030', '2040', '2050']],
        co2_df, sum_df, non_co2_df])
    summary_df['percch_2030'] = (
        summary_df['2030'] - summary_df['2020']) / summary_df['2020']
    summary_df['percch_2040'] = (
        summary_df['2040'] - summary_df['2020']) / summary_df['2020']
    summary_df['percch_2050'] = (
        summary_df['2050'] - summary_df['2020']) / summary_df['2020']
    summary_df.to_csv(
        os.path.join(_OUT_DIR, '20231117_cs_benchmark_summary.csv'),
        index=False)


def cross_sector_benchmarks_May_2024():
    """Calculate cross sector benchmarks from AR6 and IEA NZE.

    This analysis was followed for the draft of the cross-sector pathway
    revision that was finalized by Ginger in May 2024. It includes
    harmonization with historical emissions.
    
    """
    # This performs calculation of the scenario set for gross fossil
    # CO2 emissions, including harmonization with historical emissions within
    # a budget constraint calculated over 2022-2050.
    # The output is written to file at the path
    # os.path.join(_OUT_DIR, "combined_em_df_budg2050.csv")
    # summarize_CO2()

    # Gross fossil CO2 (output from above)
    harm_df = pandas.read_csv(
        os.path.join(_OUT_DIR, "combined_em_df_budg2050.csv"))

    # fill missing data via interpolation
    year_col = [col for col in harm_df if col.startswith('2')]
    harm_df.loc[
        harm_df['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross|Harmonized-DB',
        '2020'] = _2020_CO2
    harm_df['2021'] = numpy.nan
    gr_co2_df = harm_df.loc[
        (harm_df['Variable'] ==
        'Emissions|CO2|Energy and Industrial Processes|Gross|Harmonized-DB') |
        (harm_df['scen_id'] == 'NZE 2023')]
    gr_co2_df.replace(0, numpy.nan, inplace=True)
    gr_co2_df_interp = gr_co2_df[year_col].interpolate(axis=1)

    # median of harmonized scenarios, plus NZE
    gr_co2_med_ser = gr_co2_df_interp.quantile(q=0.5)
    gr_co2_med = pandas.DataFrame(gr_co2_med_ser).transpose()

    # summarize single-gas pathways for non-CO2 GHGs
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    num_cols = ['2020', '2030', '2040', '2050']
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']

    eip_n2o_df = calc_eip_n2o(ar6_scen, year_col)
    eip_n2o_df.reset_index(inplace=True)
    med_c1_n2o = eip_n2o_df.loc[
        eip_n2o_df['scen_id'].isin(c1_scen)][num_cols].quantile(q=0.5)
    n2o_df = pandas.DataFrame(med_c1_n2o).transpose()
    n2o_df['Variable'] = 'Fossil N2O'
    n2o_df['Source'] = 'Cross sector benchmark'

    eip_ch4_df = calc_eip_ch4(ar6_scen, year_col)
    eip_ch4_df.reset_index(inplace=True)
    med_c1_ch4 = eip_ch4_df.loc[
        eip_ch4_df['scen_id'].isin(c1_scen)][num_cols].quantile(q=0.5)
    ch4_df = pandas.DataFrame(med_c1_ch4).transpose()
    ch4_df['Variable'] = 'Fossil CH4'
    ch4_df['Source'] = 'Cross sector benchmark'

    med_c1_hfc = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|HFC') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)
    med_c1_pfc = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|PFC') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)
    med_c1_sf6 = ar6_scen.loc[
        (ar6_scen['Variable'] == 'Emissions|SF6') &
        (ar6_scen['scen_id'].isin(c1_scen)), ][num_cols].quantile(q=0.5)

    # add CO2eq from non-CO2 GHGs
    non_co2_df = pandas.DataFrame(
        {'N2O (Mt CO2e)': med_c1_n2o * _N2O_GWP100_AR6 * _KT_to_MT,
        'CH4 (Mt CO2e)': med_c1_ch4 * _CH4FOSS_GWP100_AR6,
        'HFC (Mt CO2e)': med_c1_hfc * _HFC_GWP100_AR6 * _KT_to_MT,
        'PFC (Mt CO2e)': med_c1_pfc * _PFC_GWP100_AR6 * _KT_to_MT,
        'SF6 (Mt CO2e)': med_c1_sf6 * _SF6_GWP100_AR6 * _KT_to_MT}).transpose()
    non_co2_df['Variable'] = non_co2_df.index
    non_co2_df['Source'] = 'Cross sector benchmark'

    co2_df = gr_co2_med[num_cols]
    co2e_df = pandas.concat([co2_df, non_co2_df])
    sum_ser = co2e_df[num_cols].sum()
    sum_df = pandas.DataFrame(sum_ser).transpose()
    sum_df['Variable'] = 'Gross fossil CO2e'
    sum_df['Source'] = 'Cross sector benchmark'

    # summary for inclusion in the report
    co2_df['Variable'] = 'Gross fossil CO2'
    co2_df['Source'] = 'Cross sector benchmark'
    summary_df = pandas.concat([co2_df, sum_df, non_co2_df])
    summary_df['percch_2030'] = (
        summary_df['2030'] - summary_df['2020']) / summary_df['2020']
    summary_df['percch_2040'] = (
        summary_df['2040'] - summary_df['2020']) / summary_df['2020']
    summary_df['percch_2050'] = (
        summary_df['2050'] - summary_df['2020']) / summary_df['2020']
    summary_df.to_csv(
        os.path.join(_OUT_DIR, '20240511_cs_benchmark_summary.csv'),
        index=False)


def summarize_n2o_ch4():
    """Summarize N2O and CH4 in C1 scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    sum_cols = year_col + ['scen_id']
    summary_cols = ['2020', '2030', '2040', '2050']

    # N2O
    n2o_var_list = [
        'Emissions|N2O|Energy', 'Emissions|N2O|Industrial Processes',
        'Emissions|N2O|Other', 'Emissions|N2O|Waste']
    sum_cols = year_col + ['scen_id']
    n2o_df = ar6_scen.loc[ar6_scen['Variable'].isin(n2o_var_list)]
    eip_n2o_df = n2o_df[sum_cols].groupby('scen_id').sum()
    eip_n2o_df.reset_index(inplace=True)
    c1_med_n2o = eip_n2o_df.loc[
        eip_n2o_df['scen_id'].isin(c1_scen)][summary_cols].quantile(q=0.5)

    # CH4
    ch4_var_list = [
        'Emissions|CH4|Energy', 'Emissions|CH4|Industrial Processes',
        'Emissions|CH4|Other', 'Emissions|CH4|Waste']
    sum_cols = year_col + ['scen_id']
    ch4_df = ar6_scen.loc[ar6_scen['Variable'].isin(ch4_var_list)]
    eip_ch4_df = ch4_df[sum_cols].groupby('scen_id').sum()
    eip_ch4_df.reset_index(inplace=True)
    c1_med_ch4 = eip_ch4_df.loc[
        eip_ch4_df['scen_id'].isin(c1_scen)][summary_cols].quantile(q=0.5)
    print('break')


def summarize_2030_renewables():
    """Summarize deployment of renewables in 2030."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)
    c1_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_scen)]
    c1_em.set_index('scen_id', inplace=True)

    filtered_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_filtered)]
    filtered_em.set_index('scen_id', inplace=True)
    summary_cols = ['2030', 'Variable']

    c1_prien_2030 = c1_em.loc[
        c1_em['Variable'] == 'Primary Energy']['2030']
    c1_renen_2030 = c1_em.loc[
        c1_em[
            'Variable'] == 'Primary Energy|Renewables (incl. Biomass)']['2030']
    c1_en_df = pandas.DataFrame({
        'Primary Energy': c1_prien_2030,
        'Primary Energy|Renewables': c1_renen_2030})
    c1_ren_share_2030 = (c1_en_df['Primary Energy|Renewables'] /
        c1_en_df['Primary Energy'])
    c1_25perc = c1_ren_share_2030.quantile(q=0.25)
    c1_75perc = c1_ren_share_2030.quantile(q=0.75)

    filt_prien_2030 = filtered_em.loc[
        filtered_em['Variable'] == 'Primary Energy']['2030']
    filt_renen_2030 = filtered_em.loc[
        filtered_em[
            'Variable'] == 'Primary Energy|Renewables (incl. Biomass)']['2030']
    filt_en_df = pandas.DataFrame({
        'Primary Energy': filt_prien_2030,
        'Primary Energy|Renewables': filt_renen_2030})
    filt_ren_share_2030 = (filt_en_df['Primary Energy|Renewables'] /
        filt_en_df['Primary Energy'])
    filt_med = filt_ren_share_2030.quantile(q=0.5)


def summarize_ip_ccs():
    """Summarize industrial process CCS from filtered scenarios."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)
    c1_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_scen)]
    c1_em.set_index('scen_id', inplace=True)

    filtered_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_filtered)]
    ip_ccs = filtered_em.loc[
        filtered_em['Variable'] ==
        'Carbon Sequestration|CCS|Industrial Processes']
    filt_med = ip_ccs['2050'].quantile(q=0.5)
    print(filt_med)


def summarize_kyoto_gases():
    """Summarize median of C1 scenarios for all Kyoto Protocol gases."""
    kyoto_var_list = [
        'Emissions|HFC',  # kt HFC134a-equiv/year
        'Emissions|PFC',  # kt CF4-equiv/year
        'Emissions|SF6']  # kt SF6/year
    ar6_key, ar6_scen = read_ar6_data()
    summary_cols = ['2020', '2030', '2040', '2050']

    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)
    c1_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_scen)]
    filtered_em = ar6_scen.loc[ar6_scen['scen_id'].isin(c1_filtered)]

    kp_c1_df = c1_em.loc[c1_em['Variable'].isin(kyoto_var_list)]
    num_scen = kp_c1_df.groupby('Variable').count()
    kp_C1_med = kp_c1_df.groupby('Variable')[summary_cols].quantile(q=0.5)
    kp_C1_med['num scen'] = num_scen['2050']
    kp_C1_med['source'] = 'median of C1'

    kp_filt_df = filtered_em.loc[filtered_em['Variable'].isin(kyoto_var_list)]
    num_filt_scen = kp_filt_df.groupby('Variable').count()
    kp_filt_med = kp_filt_df.groupby('Variable')[summary_cols].quantile(q=0.5)
    kp_filt_med['num scen'] = num_filt_scen['2050']
    kp_filt_med['source'] = 'median of filtered scenarios'

    year_col = [col for col in ar6_scen if col.startswith('2')]
    eip_n2o_df = calc_eip_n2o(ar6_scen, year_col)
    eip_n2o_df.reset_index(inplace=True)
    c1_eip_n2o = eip_n2o_df.loc[eip_n2o_df['scen_id'].isin(c1_scen)]
    med_c1_n2o = c1_eip_n2o[summary_cols].quantile(q=0.5)
    n2o_c1_df = pandas.DataFrame(med_c1_n2o).transpose()
    n2o_c1_df['Variable'] = 'EIP N2O'
    n2o_c1_df['num scen'] = c1_eip_n2o.shape[0]
    n2o_c1_df['source'] = 'median of C1'

    filt_n2o = eip_n2o_df.loc[eip_n2o_df['scen_id'].isin(c1_filtered)]
    med_filt_n2o = filt_n2o[summary_cols].quantile(q=0.5)
    n2o_filt_df = pandas.DataFrame(med_filt_n2o).transpose()
    n2o_filt_df['Variable'] = 'EIP N2O'
    n2o_filt_df['num scen'] = filt_n2o.shape[0]
    n2o_filt_df['source'] = 'median of filtered scenarios'

    eip_ch4_df = calc_eip_ch4(ar6_scen, year_col)
    eip_ch4_df.reset_index(inplace=True)
    c1_eip_ch4 = eip_ch4_df.loc[eip_ch4_df['scen_id'].isin(c1_scen)]
    med_c1_ch4 = c1_eip_ch4[summary_cols].quantile(q=0.5)
    ch4_c1_df = pandas.DataFrame(med_c1_ch4).transpose()
    ch4_c1_df['Variable'] = 'EIP CH4'
    ch4_c1_df['num scen'] = c1_eip_ch4.shape[0]
    ch4_c1_df['source'] = 'median of C1'

    filt_ch4 = eip_ch4_df.loc[eip_ch4_df['scen_id'].isin(c1_filtered)]
    med_filt_ch4 = filt_ch4[summary_cols].quantile(q=0.5)
    ch4_filt_df = pandas.DataFrame(med_filt_ch4).transpose()
    ch4_filt_df['Variable'] = 'EIP CH4'
    ch4_filt_df['num scen'] = filt_ch4.shape[0]
    ch4_filt_df['source'] = 'median of filtered scenarios'

    kp_df = pandas.concat(
        [kp_C1_med, kp_filt_med, n2o_c1_df, n2o_filt_df, ch4_c1_df,
        ch4_filt_df])
    kp_df.to_csv(
        os.path.join(_OUT_DIR, '20230908_Kyoto_gases_summary.csv'))


def compare_oecd_scenarios():
    """Compare filtered to scenarios to those in OECD 2023 report."""
    ar6_key, ar6_scen = read_ar6_data()
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # AR6 scenarios described as "Paris-consistent" by Pouille et al 2023
    oecd_path = "C:/Users/ginger.kowal/Documents/Scenario review/OECD/Table AA1.csv"
    oecd_df = pandas.read_csv(oecd_path)
    oecd_df['scen_id'] = oecd_df['model'] + ' ' + oecd_df['scenario']
    oecd_scen = oecd_df['scen_id']

    year_col = [col for col in ar6_scen if col.startswith('2')]
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)

    # calculate median gross emissions for scenarios in filtered sets
    summary_cols = ['2020', '2030', '2040', '2050', 'Variable']
    med_co2_filtered_c1 = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(c1_filtered)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_filtered_c1.reset_index(inplace=True)
    med_co2_filtered_c1['Source'] = 'filtered C1'

    med_co2_oecd = ar6_gross_em.loc[
        ar6_gross_em['scen_id'].isin(oecd_scen)][summary_cols].groupby(
            'Variable').quantile(q=0.5)
    med_co2_oecd.reset_index(inplace=True)
    med_co2_oecd['Source'] = 'OECD'
    save_df = pandas.concat([med_co2_filtered_c1, med_co2_oecd])
    save_df.to_csv(
        os.path.join(_OUT_DIR, 'gross_eip_co2_filtered_vs_oecd.csv'),
        index=False)


def id_ambitious_scenarios():
    """Compare 2050 ambition in terms of gross fossil CO2."""
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    # fill EIP emissions for all C1 scenarios in the database
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    ar6_gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # make table of net CO2, gross CO2, and CDR
    ccs_var_list = [
        'Carbon Sequestration|CCS|Biomass',
        'Carbon Sequestration|Direct Air Capture',
        'Carbon Sequestration|Enhanced Weathering']
    ccs_sum = ar6_scen.loc[
        ar6_scen['Variable'].isin(ccs_var_list)].groupby('scen_id').sum()
    ccs_sum['Variable'] = 'Carbon Sequestration|Sum'
    ccs_sum.reset_index(inplace=True)

    em_table = pandas.concat(
        [ar6_filled_em.loc[ar6_filled_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes'],
        ccs_sum, ar6_gross_em])

    # cumulative CDR and cumulative net emissions
    # em_table.replace(0, numpy.nan, inplace=True)
    # col_2050 = [str(idx) for idx in list(range(2020, 2051))]
    # c1_em = em_table.loc[em_table['scen_id'].isin(c1_scen)]
    # cum_sum_2050 = c1_em[col_2050].interpolate(axis=1).sum(axis=1)
    # # set column name
    # col_2100 = [str(idx) for idx in list(range(2020, 2100))]
    # cum_sum_2100 = c1_em[col_2100].interpolate(axis=1).sum(axis=1)
    # # set column name
    # sum_df = pandas.DataFrame(
    #   {'cum_sum_20-50': cum_sum_2050,
    #   'cum_sum_20-2100': cum_sum_2100})
    # errand_df = pandas.concat(
    #   [sum_df, c1_em[['scen_id', 'Variable']]], axis=1)
    # errand_df.to_csv(
    #   "C:/Users/ginger.kowal/Desktop/c1_cum_cdr_netem_2020-2050.csv")


    # summarize % reduction in gross fossil CO2 2020-2050
    c1_em = em_table.loc[em_table['scen_id'].isin(c1_scen)]
    c1_em['filtered'] = 0
    c1_em.loc[c1_em['scen_id'].isin(c1_filtered),'filtered'] = 1
    c1_em.to_csv(
        "C:/Users/ginger.kowal/Desktop/c1_net_gross_fossil_CO2_CDR.csv")


def demonstrate_aneris():
    """Demonstration of aneris for harmonization of historical emissions."""
    # sample data for aneris
    model, hist, driver = load_data()
    for scenario in driver.scenarios():
        driver.harmonize(scenario)
    harmonized, metadata, diagnostics = driver.harmonized_results()

    data = pandas.concat([hist, model, harmonized])
    df = data[data.Region.isin(['World'])]
    df = pandas.melt(df, id_vars=aneris.iamc_idx, value_vars=aneris.numcols(df),
                 var_name='Year', value_name='Emissions')
    df['Label'] = df['Model'] + ' ' + df['Variable']
    df.head()
    df.to_csv("C:/Users/ginger.kowal/Desktop/aneris_sample_harm.csv")

    # filtered scenarios
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # fill EIP emissions for all C1 scenarios in the database
    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    an_emissions_df = ar6_filled_em.loc[
        (ar6_filled_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes') &
        (ar6_filled_em['scen_id'].isin(c1_filtered))]
    hist_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'eip_co2_emissions_historical.csv')
    regions_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'regions_regions_sectors.csv')
    config_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'aneris_regions_sectors.yaml')

    aneris_dict = harmonize_to_historical(
        an_emissions_df, hist_path, regions_path, config_path)


def summarize_CO2():
    """Summarize gross EIP CO2 including harmonization with 2022 emissions."""
    # net EIP CO2 emissions, unharmonized
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    # table of scenario metadata for filtered scenarios
    key_filt = ar6_key.loc[ar6_key['scen_id'].isin(c1_filtered)]
    sum_cols = ['Model', 'Scenario', 'Literature Reference (if applicable)']
    sum_tab = key_filt[sum_cols]
    sum_tab.to_csv(
        os.path.join(_OUT_DIR, 'filtered_scenarios_summary_table.csv'),
        index=False)

    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    net_co2_df = ar6_filled_em.loc[
        (ar6_filled_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes') &
        (ar6_filled_em['scen_id'].isin(c1_filtered))]
    net_co2_df.replace(0, numpy.nan, inplace=True)

    # gross EIP CO2 emissions (unharmonized)
    gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    gross_co2_df = gross_em.loc[
        (gross_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes|Gross') &
        (gross_em['scen_id'].isin(c1_filtered))]
    gross_co2_df.replace(0, numpy.nan, inplace=True)

    # harmonize net CO2 emissions with historical
    hist_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'eip_co2_emissions_historical.csv')
    regions_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'regions_regions_sectors.csv')
    config_path = os.path.join(
        _PROJ_DIR, 'aneris_inputs', 'aneris_regions_sectors.yaml')
    aneris_dict = harmonize_to_historical(
        net_co2_df, hist_path, regions_path, config_path)

    # calculate gross EIP CO2 emissions from harmonized net emissions
    harm_df = aneris_dict['harmonized']
    harm_net_co2_df = harm_df.loc[
        harm_df['Variable'] == 'p|Emissions|CO2|Harmonized-DB']
    harm_net_co2_df['scen_id'] = harm_net_co2_df['Scenario']
    harm_net_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes')
    non_net_df = ar6_filled_em.loc[
        (ar6_filled_em['Variable'] !=
            'Emissions|CO2|Energy and Industrial Processes') &
        (ar6_filled_em['scen_id'].isin(harm_df['Scenario'].unique()))]
    non_net_df.replace(0, numpy.nan, inplace=True)
    comb_df = pandas.concat([harm_net_co2_df, non_net_df])
    gross_harm_co2_df = calc_gross_eip_co2(comb_df, year_col, fill_df=True)
    gross_harm_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes|Gross|Harmonized-DB')
    harm_net_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes|Harmonized-DB')

    # gross EIP CO2 from IEA Net Zero Emissions by 2050 scenario
    iea_path = os.path.join(
        _PROJ_DIR, 'IEA', 'NZE_2023', 'gross_eip_co2_emissions.csv')
    iea_co2_df = pandas.read_csv(iea_path)
    iea_co2_df['scen_id'] = iea_co2_df['Source']
    iea_co2_df['Variable'] = (
        'Emissions|CO2|Energy and Industrial Processes|Gross')

    # summarize envelope
    # add iea to unharmonized and harmonized emissions df
    em_df = pandas.concat(
        [net_co2_df, harm_net_co2_df, gross_co2_df, gross_harm_co2_df,
        iea_co2_df])
    summary_cols = [
        str(idx) for idx in list(range(2020, 2051))] + ['Variable', 'scen_id']
    em_df = em_df[summary_cols]
    em_df.to_csv(os.path.join(_OUT_DIR, "combined_em_df_budg2050.csv"))


def summarize_filtered_CO2():
    """Summarize gross EIP CO2 in filtered scenarios and all C1."""
    # net EIP CO2 emissions, unharmonized
    ar6_key, ar6_scen = read_ar6_data()
    year_col = [col for col in ar6_scen if col.startswith('2')]
    num_cols = ['2020', '2030', '2040', '2050', '2060']
    c1_scen = ar6_key.loc[ar6_key['Category'] == 'C1']['scen_id']
    c1_filtered = sustainability_filters(c1_scen, ar6_scen, filter_flag=7)

    ar6_filled_em = fill_EIP_emissions(ar6_scen, c1_scen)
    net_co2_df = ar6_filled_em.loc[
        (ar6_filled_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes') &
        (ar6_filled_em['scen_id'].isin(c1_filtered))]
    net_co2_df.replace(0, numpy.nan, inplace=True)

    # gross EIP CO2 emissions (unharmonized)
    gross_em = calc_gross_eip_co2(ar6_filled_em, year_col)
    gross_co2_df = gross_em.loc[
        (gross_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes|Gross') &
        (gross_em['scen_id'].isin(c1_scen))][num_cols + ['scen_id']]
    gross_co2_df.replace(0, numpy.nan, inplace=True)

    # median of filtered scenarios
    gross_co2_filt = gross_em.loc[
        (gross_em['Variable'] ==
            'Emissions|CO2|Energy and Industrial Processes|Gross') &
        (gross_em['scen_id'].isin(c1_filtered))]
    gross_co2_filt.replace(0, numpy.nan, inplace=True)
    gr_co2_med_ser = gross_co2_filt[num_cols].quantile(q=0.5)
    gr_co2_med = pandas.DataFrame(gr_co2_med_ser).transpose()
    gr_co2_med['Variable'] = 'Emissions|CO2|Energy and Industrial Processes|Gross'
    gr_co2_med['scen_id'] = 'Median of filtered scenarios'

    comb_df = pandas.concat([gross_co2_df, gr_co2_med])
    comb_df.to_csv(
        os.path.join(_OUT_DIR, 'gross_co2_summary.csv'), index=False)


def main():
    # cross_sector_sr15()
    # filter_AR6_scenarios()
    # extract_imps()
    # iisd_filter_variations()
    # compare_ar6_filters()
    # export_data_for_fig()
    # afforestation_test()
    # summarize_final_energy()
    # afolu_co2e_ngfs()
    # summarize_c1_key_var()
    # cross_sector_benchmarks()
    # summarize_filtered_key_var()
    # summarize_n2o_ch4()
    # summarize_2030_renewables()
    # summarize_ip_ccs()
    # summarize_kyoto_gases()
    # compare_oecd_scenarios()
    # id_ambitious_scenarios()
    # summarize_CO2()
    summarize_filtered_CO2()
    # cross_sector_benchmarks_May_2024()


if __name__ == '__main__':
    main()