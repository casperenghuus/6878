#! /usr/bin/env python

import pandas as pd
import csv

conversion = pd.read_csv('../data/ensg_to_sym.tsv', sep = '\t')
table = {}
for row in conversion.itertuples():
    table[row[1]] = row[2]
df = pd.read_csv('../data/Reference/exp_inweb_im_ensg.csv', header = None, sep = ' ')
converted = df.applymap(lambda x: table.get(x))
converted = converted.dropna()
converted = converted.rename(columns = {0:'V0', 1:'V1'})
converted.to_csv('../data/Reference/exp_inweb_im.csv', quoting = csv.QUOTE_ALL)
