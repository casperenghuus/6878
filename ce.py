#!/usr/bin/env python

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

f_Loayza = 'data/GSE42509_PQST_rna_rp_fpkms.txt'
CWD = os.getcwd()

def read_Loayza():
	'''Load data from Loayza paper into DataFrame'''
	f = os.path.join(CWD, f_Loayza)
	df = pd.read_table(f)
	return df

def heatmap(df, fn, L=None):
	'''Plot heatmap of dataframe'''
	# Make gene names into
	df = df.set_index('sym')
	# Sort by sum of rows
	idx = df.sum(axis=1).sort_values(ascending=False).index

	# Plot
	fig, ax = plt.subplots(figsize=(10,30))
	# Plot subset of data?
	if L != None:
		sns.heatmap(df.ix[idx][0:L])
	else:
		sns.heatmap(df.ix[idx])
	# Annotations
	plt.xticks(rotation=45)
	plt.yticks(rotation=0)
	plt.ylabel('Gene Symbol')
	plt.xlabel('Experiment')

	# Save figure
	plt.savefig(fn)

df = read_Loayza()
heatmap(df.drop(['enstxids'], axis=1), 'heatmap.png', 100)


