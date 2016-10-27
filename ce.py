#!/usr/bin/env python

import os
import subprocess
import pandas as pd

f_Loayza = 'data/GSE42509_PQST_rna_rp_fpkms.txt.gz'
CWD = os.getcwd()

def read_Loayza():
	# Load path and escape characters
	f = os.path.join(CWD, f_Loayza)
	# f = f.replace('(', '\\(').replace(')', '\\)').replace(' ', '\\ ')

	# Pipe .gz file
	command = 'gunzip %s'%f
	handle = subprocess.Popen(command, shell=True)
	df = pd.read_table(f[:-3])
	print df.head(50)

df = read_Loayza()


