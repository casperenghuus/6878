#!/usr/bin/env python

import os
import pandas as pd

f_Loayza = 'data/GSE42509_PQST_rna_rp_fpkms.txt'
CWD = os.getcwd()

def read_Loayza():
	'''Load data from Loayza paper into DataFrame'''
	f = os.path.join(CWD, f_Loayza)
	df = pd.read_table(f)
	return df

df = read_Loayza()


