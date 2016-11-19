#!/usr/bin/env python

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from sklearn.cluster.bicluster import SpectralBiclustering as bc
from sklearn.cluster.bicluster import SpectralCoclustering as cc
from sklearn.metrics import consensus_score

f_Loayza = 'data/GSE42509_PQST_rna_rp_fpkms.txt'
CWD = os.getcwd()

def read_Loayza():
    '''Load data from Loayza paper into DataFrame'''
    f = os.path.join(CWD, f_Loayza)
    df = pd.read_table(f)
    cols = list(df.columns.values)
    # Switch column 1 and 2
    cols_reordered = [
        cols[1],
        cols[0],
        cols[2],
        cols[8],
        cols[3],
        cols[9],
        cols[4],
        cols[10],
        cols[5],
        cols[11],
        cols[6],
        cols[12],
        cols[7],
        cols[13],
        ]
    # cols2 = [cols[1], cols[0]] + cols[2:]
    df = df[cols_reordered]

    # Add ln-transformed columns to df
    df = df.merge(
        ln_values(df.drop(['sym','enstxids'], axis=1)),
        left_index=True, right_index=True)

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

def ln_values(df):
    '''Return new columns holding ln of each element across
    all columns in df. Add 'ln' prefix to columns.'''

    # In case of negative values, take log of positive value and change sign
    func = lambda f: np.sign(f) * np.log(abs(f))

    col_names = list(df.columns.values)
    new_names = ['ln_'+name for name in col_names]
    df = df.applymap(func)
    df.columns = new_names
    return df

def Spectral_BiClustering(data, row_idx, col_idx, nClusters, prefix=''):
    '''Function to perform bipartite clustering'''
    # Create model
    model = bc(n_clusters=nClusters)
    model.fit(data)

    # Fit to data
    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]

    # Heatmap of original data
    plot_spectral(
        {'data':data,'x':col_idx,'y':row_idx},
        {'title':'Original Data','filename':prefix+'original_data.png'})

    # Heatmap of fitted data
    plot_spectral(
        {'data':fit_data,'x':col_idx,'y':row_idx},
        {'title':"After biclustering; rearranged to show biclusters",
            'filename':prefix+'bipartite_clustering.png'})

def Spectral_CoClustering(data, row_idx, col_idx, nClusters, prefix=''):
    '''Function to perform bipartite clustering'''
    # Create model
    model = cc(n_clusters=nClusters)
    model.fit(data)

    # Fit to data
    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]

    # Heatmap of fitted data
    plot_spectral(
        {'data':fit_data,'x':col_idx,'y':row_idx},
        {'title':"After co-clustering; rearranged to show co-clusters",
            'filename':prefix+'co_clustering.png'})

def plot_spectral(data, strings, log_cmap=False):
    '''Function to plot bipartite cluster'''
    plt.figure(figsize=(10,20))
    # Toggle between log-scaled colormap
    if log_cmap:
        sns.heatmap(data['data'],
            xticklabels=data['x'],
            norm=LogNorm(vmin=data.min(), vmax=data.max()))
    else:
        sns.heatmap(data['data'], xticklabels=data['x'])

    plt.xticks(rotation='vertical')
    plt.title(strings['title'])
    plt.savefig(strings['filename'])
    plt.clf(); # Clear all plots


if __name__=='__main__':
    df = read_Loayza()

    row_idx = df['sym']

    # Log-transformed columns
    prefix = 'ln_'
    col_idx = [col for col in df.columns if 'ln_' in col]
    data = df[col_idx].as_matrix()
    data[data==0] = 0.001 # Clustering does not accept zeros

    # Heatmap of original data
    plot_spectral(
        {'data':data,'x':col_idx,'y':row_idx},
        {'title':'Original Data','filename':prefix+'original_data.png'})

    # Clustering
    Spectral_BiClustering(data, row_idx, col_idx, 6, prefix)
    Spectral_CoClustering(data, row_idx, col_idx, 6, prefix)

    # Columns with raw data
    # col_idx = df.columns[2:14]
    # data = df[col_idx].as_matrix()
    # bipartite_clustering(data, row_idx, col_idx)






    # DEPRECATED COMMANDS:
    # Save Loayza data with column 0 and 1 switched
    # df.to_csv('data/GSE42509_PQST_rna_rp_fpkms_c0c1switch.txt',
    #     '\t', index=False)

    # Plot heatmap
    # heatmap(df.drop(['enstxids'], axis=1), 'heatmap.png', 100)


