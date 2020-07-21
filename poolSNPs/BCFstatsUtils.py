# Utils for plotting output of bcftools stats
import csv
import json
import pandas as pd

"""
Inspired from code of © Kristiina Ausmees
"""

###
#
# Get all rows from a file-
#
###
def get_all_rows(file):
	with open(file, 'r') as csvfile:
		rows = []
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			rows.append(row)
	return rows

###
#
# Get all rows from a file that are between a row that starts with '# id' and the next row that starts with '#'
#
# e.g. file = '/media/kristiina/My Passport/Data/Mesolithic/results/concordance/ans17/stats_HCimp_1x_only_ans17.output'
#      id = 'GCsAF'
###

def get_table(file, id):
	table_rows = []
	found = False

	with open(file, 'r') as csvfile:
		reader = csv.reader(csvfile, delimiter='\t')
		for row in reader:
			if row[0] == ("# " + id):
				found = True
				table_rows.append(row)
				continue
			if found:
				if row[0].startswith("#"):
					break
				table_rows.append(row)

	# if (id =="GCsAF"):
	# 	print("table rows")
		# print(table_rows)
	return table_rows

###
#
# Get all rows from a file that are between a row that starts with '# id' and the next row that starts with '#'
# as a dictionary
#
# Only works for id = NRDs GCTs or GCsS
#
# Dictionary is stored column-wise (1 key per colum name and list of values over all table rows. All values as str.
#
###
def get_table_dict(file, id):
	table_rows = get_table(file, id)
	print(table_rows)
	cols_to_skip = {"# NRDs": 2, "# GCTs": 2, "# GCsS": 2, "# GCsAF": 2, '# AF': 2, '# SN': 2}
	data = {}

	colid = table_rows[0][0]

	labels = [l.split("]")[1] for l in table_rows[0][cols_to_skip[colid]:]]
	for i in range(len(labels)):
		data[labels[i]] = []

	for trow in table_rows[1:]:
		values = trow[cols_to_skip[colid]:]
		map(float, values)
		# print("\n\n")
		# print(file)
		for i in range(len(labels)):
			data[labels[i]].append(values[i])
			#print(labels[i] + " : " + str(values[i]))

	return {id : data}


def get_table_dataframe(file, id):
	table_dict = get_table_dict(file, id)
	df = pd.DataFrame.from_dict(table_dict[id], orient='columns')

	return df


if __name__ == '__main__':
	### testing ###

	import os
	import matplotlib.pyplot as plt

	os.chdir('/home/camille/PoolImpHuman/data/omniexpress')
	filestats = os.path.join(os.getcwd(), 'filestats.vchk')
	gcts = get_table_dict(filestats, "GCTs")
	gcss = get_table_dataframe(filestats, "GCsS")

	gcsaf = get_table_dataframe(filestats, "GCsAF")
	af = get_table_dataframe(filestats, "AF")

	myaf = af.iloc[:, 1].to_frame().join(gcsaf)
	myaf.drop('number of genotypes', axis=1, inplace=True)
	myaf = myaf.astype(float)
	#myaf['percent SNPs'] = myaf['number of SNPs'].apply(lambda x: x / myaf['number of SNPs'].sum())

	myaf.plot()
	plt.show()