# Reproduce aneris examples
import pandas
import seaborn as sns
import matplotlib.pyplot as plt

import aneris
from aneris.tutorial import load_data

# %matplotlib inline


def getting_started():
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

	sns.lineplot(x=df.Year.astype(int), y=df.Emissions, hue=df.Label)
	plt.legend(bbox_to_anchor=(1.05, 1))


def main():
	getting_started()


if __name__ == '__main__':
    main()