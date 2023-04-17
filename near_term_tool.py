"""Replicate functionality of near term target setting tool."""
import os
import tempfile
import shutil

import numpy
import pandas

# back end data dictionary
_BASE_DATA_PATH_DICT = {
    # placeholder
    # TODO: decide how you want to store and load
                              # backend data like sector activity and emissions.
                              # I recommend you store them in external csv
                              # files stored in the repository, like we did for
                              # buildings
    'power': {
    	'activity': '<path to csv file>',
    	'emissions': '<path to csv file>',
    },
    # add other sectors
}


# available target setting methods
_METHOD_LIST = ['SDA', 'ARA']

# available scenarios
_SCENARIO = 'SBTi_1.5C'

# available sectors
_SDA_SECTOR_LIST = ['power', 'buildings_res', 'buildings_serv', 'cement']

# available activity types
_ACTIVITY_TYPES = ['fixed_market_share', 'ty_output']

# directory to hold outputs
outdir = 'dir'  # placeholder


def main():
    """Execute workflow."""
    args_dict = construct_args_dict()
    sda_target_dict = near_term_target('SDA', args_dict)
    ara_target_dict = near_term_target('ARA', args_dict)


def construct_args_dict():
    """Construct a dictionary of args to be used as input to tool."""
    # This currently returns a single dictionary of default inputs.
    # As you expand the scope of your sensitivity analysis, you can modify
    # the inputs to this function to return a dictionary of inputs with values
    # that have been changed according to some sampling scheme or sensitivity
    # analysis design.
    args = {
        'scenario' = _SCENARIO;
        'sector': 'power',
        'base_year': 2020,
        'by_activity': 100,
        'by_scope1': 100,
        'by_scope2': 'NA',
        'target_year': 2030,
        'ty_activity_type': 'fixed_market_share',
        'ty_activity': 'NA',
        'most_recent_year': 'NA',
        'mry_scope1': 'NA',
        'mry_scope2': 'NA',
    }
    return args


def calc_SDA_target(args_dict):
    """
    Calculate near-term target using the SDA.

    <description of how this function works>

    Args:
        args_dict (dictionary): dictionary containing key:value pairs for
            the following required keys:
                sector
                base_year: base year
                by_activity: base year activity
                by_scope1: base year scope 1 emissions
                by_scope2: base year scope 2 emissions
                target_year: target year
                ty_activity_type: target year activity type
                ty_activity: target year activity

    Returns:
        a dictionary containing key:value pairs for the following keys:
            ty_scope1: target year scope 1 emissions
            ty_scope2: target year scope 2 emissions

    """
    # recommended: use error checking to make sure that the values for
    # args_dict['sector'] and args_dict['ty_activity_tye'] match a valid value
    # in _SDA_SECTOR_LIST and _ACTIVITY_TYPES respectively

    # use args_dict['sector'] to retrieve the right sector backend data
    sector_activity = _BASE_DATA_PATH_DICT[args_dict['sector']]['activity']
    sector_emissions = _BASE_DATA_PATH_DICT[args_dict['sector']]['emissions']

    # use the backend data and values in args_dict to calculate SDA target
    # you can copy and then modify the code from our buildings script in R!
    # https://github.com/vakowal/sda-regional-pathways/blob/main/scripts/00_packages_functions_globals.R

    # return a dictionary that looks like this:
    # target_dict = {
    #   'ty_scope1': target_year_scope1,
    #   'ty_scope2': target_year_scope2,
    #}
    # return target_dict
    pass  # placeholder


def calc_ARA_target(args_dict):
    """
    Calculate near-term target using the Absolute Reduction Approach.

    <description>

	Args:
        args_dict (dictionary): dictionary containing key:value pairs for
            the following required keys:
                base_year: base year
                by_scope1: base year scope 1 emissions
                by_scope2: base year scope 2 emissions
                target_year: target year
                most_recent_year: most recent year
                mry_scope1: scope 1 emissions in most recent year
                mry_scope2: scope 2 emissions in most recent year

    Returns:
        a dictionary containing key:value pairs for the following keys:
            ty_scope1: target year scope 1 emissions
            ty_scope2: target year scope 2 emissions

    """
    # use the values in args_dict to calculate ARA target

    # return a dictionary that looks like this:
    # target_dict = {
    #   'ty_scope1': target_year_scope1,
    #   'ty_scope2': target_year_scope2,
    #}
    # return target_dict
    pass  # placeholder


def near_term_target(method, args_dict):
    """Calculate near term target.

    <description of the function>

    Args:
        method (string): <parameter description>
        <param> (<param type>): <description>

    Returns:
        <value returned by the function>

    """
    if method == 'SDA':
        calc_SDA_target()

    elif method == 'ARA':
        calc_ARA_target()

    else:
        raise ValueError()  # TODO: best practice is to include an informative
                            # message here about the error


if __name__ == '__main__':
    main()
