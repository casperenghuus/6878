#!/usr/bin/env python

import os
import re
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import BipartiteClustering as bpc


# Path to parent directory of tp files
pipelinePath = 'FortunatoPipelinePath.txt'
temp_path = 'data/Examples/tp_enr.csv'

def build_df(args):
    '''Build matrix from tp column vectors'''

    df = None

    for s in args.slices:
        try:
            for p in args.paths[s]:
                if p['cluster'] == '1':
                    modules = bpc.tp_to_col(p['path'], args)

        except KeyError, e:
            print '\nWarning: %s not a valid slice\n'%e

    # Convert dict to dataframe
    df = pd.DataFrame(modules).transpose()
    df['module'] = modules.keys()

    cols = df.drop('nodes', axis=1).columns.values
    df[cols] = df[cols].apply(pd.to_numeric)
    return df

def build_enr(args):
    '''Build df from enrichment file'''

    df = pd.read_csv(args.enr_path)
    df.rename(columns={'Unnamed: 0':'module'}, inplace=True)

    chr_cols = [c for c in df.columns if 'chr' in c]
    pthw_cols = [c for c in df.columns if 'BIOCARTA' in c]

    # Number of features present / number of significant features present
    df['n_chr_features'] = count_feats(df[chr_cols], 1.0)
    df['n_sign_chr_features'] = count_feats(df[chr_cols], args.pthres, True)

    df['n_pthw_features'] = count_feats(df[pthw_cols], 1.0)
    df['n_sign_pthw_features'] = count_feats(df[pthw_cols], args.pthres, True)

    # Merge with existing df and return
    df = args.df.merge(df, on='module')
    return df

def count_feats(df, thres, le=False):
    '''Return number of features below given threshold'''
    if le:
        return df[df <= thres].count(axis=1)
    else:
        return df[df < thres].count(axis=1)

def plot_enr(args):
    '''Plot enrichment results'''
    df = args.df.copy()
    dfc = df[(df['n_chr_features'] > 0) & (df['size'] > 1)]
    dfp = df[(df['n_pthw_features'] > 0) & (df['size'] > 1)]

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        nrows=2, ncols=2, figsize=args.psize)

    # PLOT 1
    ax1.scatter(
        x=dfc['size'],
        y=dfc['n_chr_features'],
        label='Chromosome features',
        c='k')
    ax1.scatter(
        x=dfc['size'],
        y=dfc['n_sign_chr_features'],
        label='Significant chromosome features',
        c='r')

    ax1.set_xlabel('Number of genes in module')
    ax1.set_ylabel('Number of chromosome features in module')
    ax1.legend()

    # PLOT 2
    ax2.scatter(
        x=dfc['module'],
        y=dfc['n_sign_chr_features']/dfc['n_chr_features'],
        label='Sized by # of members in module',
        s=dfc['size'])
    ax2.set_xlabel('Module')
    ax2.set_ylabel('n Sign. Chr. features / n Chr. features')
    ax2.legend()

    # PLOT 3
    ax3.scatter(
        x=dfp['size'],
        y=dfp['n_pthw_features'],
        label='Pathway features',
        c='k')
    ax3.scatter(
        x=dfp['size'],
        y=dfp['n_sign_pthw_features'],
        label='Significant pathway features',
        c='r')

    ax3.set_xlabel('Number of genes in module')
    ax3.set_ylabel('Number of pathway features in module')
    ax3.legend()

    # PLOT 4
    ax4.scatter(
        x=dfp['module'],
        y=dfp['n_sign_pthw_features']/dfp['n_pthw_features'],
        label='Sized by # of members in module',
        s=dfp['size'])
    ax4.set_xlabel('Module')
    ax4.set_ylabel('n Sign. Pathw. features / n Pathw. features')
    ax4.legend()

    plt.savefig(args.pf + args.pname)

def plot_promiscuity(df, prefix, args):
    df['promiscuity'] = count_feats(df, 1.0)
    df['sign_promiscuity'] = count_feats(df, args.pthres)
    df = df.sort_values(
        by=['promiscuity', 'sign_promiscuity'],
        ascending=[False, False])

    if prefix == 'pthw_':
        labels = []
        for l in df.index.values:
            labels.append(l.split('_')[1])

    fig = plt.figure(figsize=(args.psize[0], args.psize[1]/4))

    sns.barplot(
        x=df.index,
        y=df['promiscuity'],
        label='Feature is present',
        color='k')
    if prefix == 'pthw_':
        sns.barplot(
            x=df.index,
            y=df['sign_promiscuity'],
            label='Feature is present and significant',
            color='r')
        plt.set(xticklabels=labels)
    else:
        sns.barplot(
            x=df.index,
            y=df['sign_promiscuity'],
            label='Feature is present and significant',
            color='r')

    plt.legend()
    plt.xticks(rotation='vertical', size='x-small')
    plt.ylabel('number of modules')
    if prefix == 'chr_':
        plt.xlabel('Chromosome')
    else:
        plt.xlabel('Pathway')

    plt.savefig(args.pf + prefix + args.pname2)
    plt.clf

if __name__=='__main__':

    # Get file paths and number of clusters
    paths, nClusters = bpc.build_paths()

    # PARSE ARGUMENTS
    parser = argparse.ArgumentParser(description = 'Parse input network files')
    parser.add_argument(
        '-s',
        dest = 'slices',
        nargs = '+',
        help = ' '.join([
            'Slices to work on.',
            'Takes space-separated inputs, i.e. "1 3 14".',
            '(Default = all)']),
        default = ['All'],
        type = str)
    parser.add_argument(
        '-pthres',
        dest = 'pthres',
        help = 'Threshold for declaring significance. (Default = 0.05)',
        default = 0.05,
        type = float)
    parser.add_argument(
        '-p',
        dest = 'plot',
        help = 'Set to plot heatmaps',
        nargs = '?',
        default = False,
        const = True,
        type = bool)
    parser.add_argument(
        '-pf',
        dest = 'pf',
        help = 'Plot folder to save plots. (Default = "plots/")',
        default = 'plots/',
        type = str)
    parser.add_argument(
        '-pname1',
        dest = 'pname',
        help = ' '.join([
            'Filename for output plot',
            '(Default = "tp_features")']),
        default = 'tp_features',
        type = str)
    parser.add_argument(
        '-pname2',
        dest = 'pname2',
        help = ' '.join([
            'Filename for output plot of promiscuity',
            '(Default = "tp_prmscty_features")']),
        default = 'tp_features',
        type = str)
    parser.add_argument(
        '-ps',
        dest = 'psize',
        nargs = 2,
        help = ' '.join([
            'Width, height of plot (space separated values).',
            '(Default = "20 20")']),
        default = [20, 30],
        type = int)

    args = parser.parse_args()

    # Update args
    args.paths = paths
    args.enr_path = temp_path # TEMPORARY PATH!!!
    args.set_module = True
    if 'All' in args.slices:
        args.slices = args.paths.keys()

    # BUILD DF
    args.df = build_df(args)
    args.df = build_enr(args)

    # Plot
    if args.plot:
        plot_enr(args)

        chr_cols = [c for c in args.df.columns if 'chr' in c]
        pthw_cols = [c for c in args.df.columns if 'BIOCARTA' in c]

        plot_promiscuity(args.df[chr_cols].T, 'chr_', args)
        plot_promiscuity(args.df[pthw_cols].T, 'pthw_', args)







