#! /usr/bin/env python2
import subprocess
from joblib import Parallel, delayed

def execCommand(x):
    subprocess.call(x, shell = True)

execCommand('python2 BipartiteClustering.py -fout biCl_red-new -p -po red-new_noClust_def-Clusters -pc red-new_def-Clusters')

commandlist = []

sizes = [190, 170, 150, 110, 90, 70, 50]
for size in sizes:
    commandlist.append('python2 BipartiteClustering.py -fout biCl_red-new_c{} -c 190'.format(size, size))


Parallel(n_jobs = 24)(delayed(execCommand)(x) for x in commandlist)
