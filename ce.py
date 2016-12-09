#!/usr/bin/env python

import os
import re
import random
import pygraphviz
from networkx.drawing.nx_agraph import graphviz_layout
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import seaborn as sns
from scipy.stats import linregress
from sklearn.cluster.bicluster import SpectralBiclustering as bc
from sklearn.cluster.bicluster import SpectralCoclustering as cc
from sklearn.metrics import consensus_score
import mygene

CWD = os.getcwd()
plot_folder = 'plots/'
f_Loayza = 'data/GSE42509_PQST_rna_rp_fpkms.txt'
f_Astrocytes = 'data/GSE58910_SenescentAstrocytesDEanalysis_v2.csv'
f_Arachne = [
    'data/net_tf_final.txt',
    'data/net_mirna_final.txt',
    'data/net_exp_final.txt',
    'data/nodes.txt'
    ]
f_jenage = [
    'data/GSE63577_counts_rpkm_exvivo_jenage_data.txt',
    'data/GSE63577_exVivo_counts_RPKMs_II.txt'
    ]
refseq_to_ensg = 'data/uniq_geneid_sym.tsv'
ensg_to_sym = 'data/ensg_to_sym.tsv'


def read_Loayza():
    '''Load data from Loayza paper into DataFrame'''
    f = os.path.join(CWD, f_Loayza)
    df = pd.read_table(f)
    cols = list(df.columns.values)
    # Switch column 1 and 2 and reorder by experimental condition
    # cols_reordered = [
    #     cols[1],
    #     cols[0],
    #     cols[2],
    #     cols[8],
    #     cols[3],
    #     cols[9],
    #     cols[4],
    #     cols[10],
    #     cols[5],
    #     cols[11],
    #     cols[6],
    #     cols[12],
    #     cols[7],
    #     cols[13],
    #     ]
    # Switch column 1 and 2 and reorder by experimental condition
    cols_reordered = [cols[1], cols[0]] + cols[2:]
    df = df[cols_reordered]

    # Add ln-transformed columns to df
    # df = df.merge(
    #     ln_values(df.drop(['sym','enstxids'], axis=1)),
    #     left_index=True, right_index=True)

    return df

def read_Astrocytes():
    f = os.path.join(CWD, f_Astrocytes)
    df = pd.read_csv(f)
    df_ids = pd.read_table(refseq_to_ensg, sep='\t')
    df_syms = pd.read_table(ensg_to_sym, sep='\t')
    df_syms.columns = ['gene_id', 'symJ']

    # Drop columns that do not contain gene information
    ids_to_drop = [
        'alignment_not_unique',
        'ambiguous','no_feature',
        'not_aligned',
        'too_low_aQual']
    df = df[~df['id'].isin(ids_to_drop)]

    df = pd.merge(
        df.drop(['Unnamed: 0'], axis=1),
        df_ids,
        on='id',
        how='left').merge(
            df_syms,
            on='gene_id',
            how='left')

    # Discard columns with RPKM < 1
    df = df.loc[(df[['baseMeanA', 'baseMeanB']].sum(axis=1) > 1)]

    # # If duplicate genes, i.e. from different splice
    # #  variants, merge their counts by summing
    # df[['astro_AB_Mean', 'astro_A', 'astro_B']] = df.groupby(
    #     ['symJ'])[
    #         ['baseMean', 'baseMeanA', 'baseMeanB']].transform('sum')
    # df = df[df.duplicated(['symJ'], keep=False)].sort_values(
    #     ['symJ', 'baseMean'])

    # # Drop duplicate rows to prevent duplicate gene_ids when merging with
    # #  other dfs. Keep splice variant with highest number of reads only.
    # df = df.drop_duplicates('symJ', keep='last')

    # # Take ln values of columns and discard those that are no longer
    # #  true, i.e. p-vals and baseMeans for only one splice variant
    df = df[['symJ', 'gene_id', 'baseMeanA', 'baseMeanB']]

    # Add ln-transformed columns to df
    # df = df.merge(
    #     ln_values(df.drop(['gene_id'], axis=1)),
    #     left_index=True, right_index=True)
    df = drop_duplicates(df, 'symJ')
    col_idx = [col for col in df.columns if '_sum' in col]

    return df[['symJ', 'gene_id'] + col_idx]

def read_hoare(load, fetch_sym):
    '''
    Merge all files for dataset into single dataframe. Add column with gene
    symbols identified from ensemble gene IDs. Function either loads from
    separate files and create a merged file, or load the merged file.
    '''

    # Output file
    fn = 'data/GSE72409.csv'
    # File holding gene symbols
    gene_file = 'data/GSE72409/geneid_to_sym.csv'
    # File holding gene lengths
    lengths_file = 'data/GSE72409/genes_len.tsv'

    if load:
        df1 = pd.read_csv(fn)
    else:
        # Load files
        to_load = {}
        target = re.compile('GSM\d+_(.+).txt')
        for root, dirs, files in os.walk(CWD, followlinks=True):
            for f in files:
                match = target.search(f)
                if match:
                    to_load[os.path.join(root, f)] = 'expr_'+match.group(1)

        # Create single, merged pandas dataframe
        initiate = True
        for k,v in to_load.items():
            if initiate:
                df1 = pd.read_table(k, header=None, names=['gene_id', v])
                initiate = False
            else:
                df2 = pd.read_table(k, header=None, names=['gene_id', v])
                df1 = pd.merge(df1, df2, on='gene_id', how='outer')

        # Add gene symbols
        if fetch_sym:
            mg = mygene.MyGeneInfo()
            genes = mg.getgenes(df1['gene_id'], 'symbol')
            gene_df = pd.DataFrame(genes)[['query', 'symbol']]
            gene_df.columns = ['gene_id', 'sym']

            # Save symbols to avoid loading them multiple times
            gene_df.to_csv(gene_file, index=False)
        else:
            # Load saved symbols if already fetched
            gene_df = pd.read_csv(gene_file)

        # Merge symbols with data and save to file
        df1 = pd.merge(gene_df, df1, on='gene_id', how='outer')
        df1.to_csv(fn, index=False)

    # Drop non-gene entries
    ids_to_drop = [
        'alignment_not_unique',
        'ambiguous','no_feature',
        'not_aligned',
        'too_low_aQual']
    df1 = df1[~df1['gene_id'].isin(ids_to_drop)]

    # Calculate RPKM
    # Merge transcript lengths into df
    lengths = pd.read_table(lengths_file, sep='\t')
    df1 = df1.merge(lengths, on='gene_id', how='outer')

    # Columns with RNA-seq counts
    cols = df1.columns.values[2:]

    # Calculate RPKM and update count columns
    func = lambda row: 10**9*row/(row['median_mRNA_len'])
    df1[cols] = df1[cols].apply(func, axis=1).div(df1[cols].sum(axis=0))

    df1 = df1.drop('median_mRNA_len', axis=1).merge(
        lengths,
        on='gene_id', how='outer')

    # Add ln-transformed columns to df
    # df1 = df1.merge(
    #     ln_values(df1.drop(['gene_id','sym','median_mRNA_len'], axis=1)),
    #     left_index=True, right_index=True)

    return df1.drop('median_mRNA_len', axis=1)

def read_jenage():
    '''Read JenAge data'''

    f1 = os.path.join(CWD, f_jenage[0])
    f2 = os.path.join(CWD, f_jenage[1])
    df1 = pd.read_table(f1, sep='\t')
    df2 = pd.read_table(f2, sep='\t')

    # Columns to merge the two dataframes on
    cols = [
        'ensembl_gene_id',
        'external_gene_id',
        'description',
        'gene_biotype'
        ]
    df = pd.merge(df1, df2, on=cols, how='outer')

    # Add ln-transformed columns to df
    # df = df.merge(
    #     ln_values(df.drop(cols, axis=1)),
    #     left_index=True, right_index=True)

    # Rename columns
    col_dict = {'external_gene_id': 'symJ', 'ensembl_gene_id': 'gene_id'}
    df.columns = [col_dict.get(x, x) for x in df.columns]

    return df

def read_arachne(single_df=False):
    '''Read Arachne results into df'''

    # Load file names
    f1 = os.path.join(CWD, f_Arachne[0])
    f2 = os.path.join(CWD, f_Arachne[1])
    f3 = os.path.join(CWD, f_Arachne[2])
    f4 = os.path.join(CWD, f_Arachne[3])

    # Column names
    colnames = ['node1', 'node2', 'edge']

    # Read to dataframe
    df1 = pd.read_table(f1, sep=' ', header=None, names=colnames)
    df2 = pd.read_table(f2, sep=' ', header=None, names=colnames)
    df3 = pd.read_table(f3, sep=' ', header=None, names=colnames)
    df4 = pd.read_table(f4, header=None, names=['sym','a','b'])
    df4 = df4.drop(['a','b'], axis=1)

    # Add column to hold node-symbol index. No offset
    # df4['id'] = df4.index
    # Add column to hold node-symbol index. Offset by 1
    df4['id'] = df4.index + 1

    df = {'tf':df1, 'mirna':df2, 'exp':df3, 'nodes':df4}

    # Return dict of dataframes
    return df

def arachne_clusters(df):
    '''Identify clusters by edges'''

    # If values are 1 indexed and not 0 indexed
    offset = 1

    # Convert dataframe to list of tuples holding node pairs
    subset = df[['node1', 'node2']]
    tuples = [tuple(x) for x in subset.values]

    # Identify clusters by edges
    G = nx.Graph()
    G.add_edges_from(tuples)
    clusters = []
    for connected_component in nx.connected_components(G):
        clusters.append(connected_component)

    # Number of unique nodes
    uniq_nodes = pd.concat([dfAra['exp'].node1, dfAra['exp'].node2])
    n_nodes = uniq_nodes.unique().shape[0]
    max_node = uniq_nodes.max() + offset

    # Matrix to hold clusters
    M_rows = np.zeros([max_node, len(clusters)], dtype=bool)
    # M_cols = np.zeros([max_node, len(clusters)], dtype=bool)
    for i in xrange(len(clusters)):
        for j in clusters[i]:
            M_rows[i,j - offset] = 1

    # # Plot graph colored by subgraphs
    # plt.figure(figsize=(10,10))
    # pos = graphviz_layout(G, prog="neato")
    # subG = nx.connected_component_subgraphs(G)
    # for g in subG:
    #     # Random color. List of length = # of nodes
    #     c = [random.random()] * nx.number_of_nodes(g)
    #     nx.draw(
    #         g,
    #         pos,
    #         node_size=30,
    #         node_color=c,
    #         cmap=plt.cm.Spectral,
    #         # vmin,max used to set colormap
    #         vmin=0.0,
    #         vmax=1.0)
    # plt.savefig(plot_folder + 'Arachne.png')
    # plt.clf(); # Clear all plots

    # Return list of clusters. Each cluster is a set
    return M

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

def drop_duplicates(df, group):
    '''Remove duplicate ids from Astrocyte data'''

    col_idx = df._get_numeric_data().columns

    df[[x + '_sum' for x in col_idx]] = df.groupby(
        group)[col_idx].transform('sum')

    df = df[
        df.duplicated([group], keep=False)
        ].sort_values([group])

    # Drop duplicate rows to prevent duplicate gene_ids when merging with
    #  other dfs. Keep splice variant with highest number of reads only.
    df = df.drop_duplicates(group, keep='last')

    return df

def Spectral_BiClustering(data, row_idx, col_idx, nClusters):
    '''Function to perform bipartite clustering'''
    # Create model
    model = bc(n_clusters=nClusters)
    model.fit(data)

    # Fit to data
    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]

    sorted_row_idx = [a for (b,a) in sorted(
        zip(np.argsort(model.row_labels_),row_idx))]
    sorted_col_idx = [a for (b,a) in sorted(
        zip(np.argsort(model.column_labels_),col_idx))]

    return model, fit_data, [sorted_row_idx, sorted_col_idx]

def Spectral_CoClustering(data, row_idx, col_idx, nClusters):
    '''Function to perform bipartite clustering'''
    # Create model
    model = cc(n_clusters=nClusters)
    model.fit(data)

    # Fit to data
    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]
    sorted_row_idx = [a for (b,a) in sorted(
        zip(np.argsort(model.row_labels_),row_idx))]
    sorted_col_idx = [a for (b,a) in sorted(
        zip(np.argsort(model.column_labels_),col_idx))]

    return model, fit_data, [sorted_row_idx, sorted_col_idx]

def jaccard_consensus(a,b):
    '''Get similarity between individual biclusters,
    a and b, calculated using the Jaccard index'''
    score = consensus_score(a,b)
    print("consensus score: {:.3f}".format(score))

def heatmap(df, rownames, fn, L=None):
    '''Plot heatmap of dataframe'''
    # Make gene names into row names
    df = df.set_index(rownames)
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
    plt.savefig(plot_folder + fn)

def plot_spectral(data, strings, figattr, log_cmap=False):
    '''Function to plot bipartite cluster'''
    plt.figure(figsize=figattr['size'])
    # Toggle between log-scaled colormap
    if log_cmap:
        sns.heatmap(data['data'],
            xticklabels=data['x'],
            norm=LogNorm(vmin=data.min(), vmax=data.max()))
    else:
        sns.heatmap(data['data'], xticklabels=data['x'])

    plt.xticks(rotation='vertical')
    plt.title(strings['title'])
    plt.savefig(plot_folder + strings['filename'])
    plt.clf(); # Clear all plots

def xy_linregress_plot(df,x,y,title,fn):
    '''X,Y plot'''

    # Ignore missing values and inf
    df = df.replace([np.inf, -np.inf], np.nan).dropna()
    slope, intercept, r_value, p_value, std_err = linregress(df[x], df[y])

    # Figure
    fig, ax = plt.subplots(figsize=(10,10))

    # Scatterplot with linear regression
    plt.scatter(df[x], df[y])
    plt.plot(df.index, df.index*slope + intercept)

    # Change if R^2 is plotted in upper/lower left corner
    if slope > 0:
        text_y = df[[x,y]].max().max()
    else:
        text_y = df[[x,y]].min().min()
    # Plot R^2
    plt.text(
        min(df.index)+1000,
        text_y*0.90,
        'y = {:.{prec}f}*{:.{prec}}+{:.{prec}f}\nR^2 = {:.{prec}f}'.format(
            slope, 'x', intercept, r_value, prec=3),
        fontsize=12)

    # Plot layout
    max_val = max(df[[x]].max().max(), df[[y]].max().max())
    interval = round(max_val/5, -len(str(int(max_val))) + 2)
    major_ticks = np.arange(0, max_val, interval)
    ax.set_xticks(major_ticks)
    ax.set_yticks(major_ticks)
    plt.xlim([-100,1.05*df[[x]].max().max()])
    plt.ylim([-100,1.05*df[[y]].max().max()])

    # Annotations
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.ylabel(x)
    plt.xlabel(y)
    plt.title(title)

    # Save figure
    plt.savefig(plot_folder + fn)


if __name__=='__main__':

    # Toggle loading data
    loayza = False
    hoare = True
    jenage = False
    astrocytes = False

    arachne = False
    tp = False

    # Toggle merged data
    merge_all = False
    hoare_only = True

    # Toggle analyses
    L_analysis = False
    HL_analysis = False
    AL_analysis = False


    if jenage:
        dfJ = read_jenage()

    if hoare:
        load = False
        fetch_sym = False
        dfH = read_hoare(load, fetch_sym)

    if arachne:
        # Returns dict of dataframes with values 'tf', 'mirna' and 'exp'
        dfAra = read_arachne()
        dfAra_m = arachne_clusters(dfAra['exp'])

    if astrocytes:
        dfA = read_Astrocytes()
        # Comparison of each sample
        # xy_linregress_plot(
        #    dfA[['astro_A', 'astro_B']],
        #    'astro_A',
        #    'astro_B',
        #    'Astrocytes',
        #    'astrocytes_xy.png')

    if loayza:
        dfL = read_Loayza()

        if L_analysis:
            # Log-transformed columns
            prefix = 'ln_'

            row_idx = dfL['sym']
            col_idx = [col for col in dfL.columns if 'ln_' in col]

            data = dfL[col_idx].as_matrix()
            data[data==0] = 0.001 # Clustering does not accept zeros

            # Comparison of control experiments
            xy_linregress_plot(
                dfL[['C12.rna', 'C3.rna']],
                'C12.rna',
                'C3.rna',
                'Loyaza C12 vs. C3 RNA-Seq',
                'loayza_C12,3_rna.png')
            xy_linregress_plot(
                dfL[['C12.rp', 'C3.rp']],
                'C12.rp',
                'C3.rp',
                'Loyaza C12 vs. C3 Ribo-Seq',
                'loayza_C12,3_rp.png')

            # Bi- and Co-clustering
            n_clusters = len(col_idx)
            n_clusters2 = 97
            bc_model, bc_fit, bc_idx = Spectral_BiClustering(
                data, row_idx, col_idx, n_clusters)
            cc_model, cc_fit, cc_idx = Spectral_CoClustering(
                data, row_idx, col_idx, n_clusters2)

            # HEATMAPS:
            # Heatmap of original data
            plot_spectral(
                {'data':data,'x':col_idx,'y':row_idx},
                {'title':'L: Original Data',
                    'filename':prefix+'L_original_data.png'},
                {'size':(10,20)})

            # Heatmap of biclustered data
            plot_spectral(
                {'data':bc_fit,'x':bc_idx[1],'y':bc_idx[0]},
                {'title':"L: After biclustering; rearranged to show biclusters",
                    'filename':prefix+'L_bipartite_clustering.png'},
                {'size':(10,20)})

            # Heatmap of coclustered data
            plot_spectral(
                {'data':cc_fit,'x':cc_idx[1],'y':cc_idx[0]},
                {'title':"L: After coclustering; rearranged to show coclusters",
                    'filename':prefix+'L_co_clustering.png'},
                {'size':(10,20)})

            # Compare bi-/co-clustered networks by Jaccard index
            jaccard_consensus(bc_model.biclusters_, cc_model.biclusters_)

    # MERGE ALL DATASETS INTO SINGLE DATAFRAME
    if merge_all:
        dfAll = pd.merge(dfH, dfL, on='sym', how='outer')
        dfAll = dfAll.merge(
            dfJ, on='gene_id', how='outer').merge(
                dfA, on='symJ', how='outer')

        print 'data size before filtering =', dfAll.shape, '\n'


        # LOG-TRANSFORM DATA
        dfAll = dfAll.merge(
            ln_values(dfAll._get_numeric_data()),
            left_index=True,
            right_index=True)


        # PREPARE HEATMAPS
        col_idx = [col for col in dfAll.columns if 'ln_' in col]

        # Data1 = weak subset. Keep rows with rowSum > 0
        data1 = dfAll.loc[(dfAll[col_idx].sum(axis=1) > 0)]
        # Data2 = strong subset. Keep only rows with all values > 0
        data2 = dfAll.replace([np.inf, -np.inf], np.nan).dropna()
        data2 = data2[(data2[col_idx] != 0).all(1)]

        row1_idx = data1['sym']
        row2_idx = data2['sym']
        print '\ndata1 size =', data1.shape
        print 'data2 size =', data2.shape, '\n'


        # SAVE DATAFRAMES TO FILE
        cols_to_save = ['sym', 'gene_id_x'] + col_idx
        data1[cols_to_save].replace(np.nan, 0).to_csv(
            'data/dfAll1.tsv', sep='\t', index=False)
        data2[cols_to_save].replace(np.nan, 0).to_csv(
            'data/dfAll2.tsv', sep='\t', index=False)


        # PLOT HEATMAPS OF DATAFRAMES
        plot_spectral(
            {'data':data1[col_idx].replace(np.nan, 0).as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data1',
                'filename':'dfAll1.png'},
            {'size':(20,20)})

        plot_spectral(
            {'data':data2[col_idx].replace(np.nan, 0).as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data2',
                'filename':'dfAll2.png'},
            {'size':(20,20)})


        # CALCULATE SD BY ROW
        # Subset by SD on data2. Cutoff = 1
        sd1 = data1[col_idx].std(axis=1)
        data1['sd'] = sd1.values
        plt.figure(figsize=(10,10))
        sns.distplot(sd1.values[~np.isnan(sd1.values)])
        plt.savefig('plots/dfAll1_sdDistr.png')

        # Subset by SD on data2. Cutoff = 0.5
        sd2 = data2[col_idx].std(axis=1)
        data2['sd'] = sd2.values
        plt.figure(figsize=(10,10))
        sns.distplot(sd2.values[~np.isnan(sd2.values)])
        plt.savefig('plots/dfAll2_sdDistr.png')

        print '\nShape after subsetting data1 on SD',data1[data1['sd']>1].shape
        print 'Shape after subsetting data2 on SD',data2[data2['sd']>0.5].shape


        # PLOT HEATMAPS SHOWING DATA AFTER SD-SUBSETTING
        # Subset by SD on data1. Cutoff = 1.0
        plot_spectral(
            {'data':data1[data1['sd'] > 0.5][col_idx].as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data1, SD > 0.5',
                'filename':'dfAll1_sdSubset.png'},
            {'size':(20,20)})

        plot_spectral(
            {'data':data2[data2['sd'] > 0.5][col_idx].as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data2, SD > 0.5',
                'filename':'dfAll2_sdSubset.png'},
            {'size':(20,20)})

        # SAVE DATAFRAMES TO FILE
        data1[cols_to_save][data1['sd']>1].replace(np.nan, 0).to_csv(
            'data/dfAll1_sdSubset.tsv', sep='\t', index=False)
        data2[cols_to_save][data2['sd']>0.5].replace(np.nan, 0).to_csv(
            'data/dfAll2_sdSubset.tsv', sep='\t', index=False)

    if hoare_only:
        print 'data size before filtering =', dfH.shape, '\n'

        # LOG-TRANSFORM DATA
        dfH = dfH.merge(
            ln_values(dfH._get_numeric_data()),
            left_index=True,
            right_index=True)

        # PREPARE HEATMAPS
        col_idx = [col for col in dfH.columns if 'ln_' in col]

        # Data1 = weak subset. Keep rows with rowSum > 0
        data1 = dfH.loc[(dfH[col_idx].sum(axis=1) > 0)]
        # Data2 = strong subset. Keep only rows with all values > 0
        data2 = dfH.replace([np.inf, -np.inf], np.nan).dropna()
        data2 = data2[(data2[col_idx] != 0).all(1)]

        row1_idx = data1['sym']
        row2_idx = data2['sym']
        print '\ndata1 size =', data1.shape
        print 'data2 size =', data2.shape, '\n'


        # SAVE DATAFRAMES TO FILE
        cols_to_save = ['sym', 'gene_id'] + col_idx
        data1[cols_to_save].replace(np.nan, 0).to_csv(
            'data/dfH1.tsv', sep='\t', index=False)
        data2[cols_to_save].replace(np.nan, 0).to_csv(
            'data/dfH2.tsv', sep='\t', index=False)


        # PLOT HEATMAPS OF DATAFRAMES
        plot_spectral(
            {'data':data1[col_idx].replace(np.nan, 0).as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data1',
                'filename':'dfH1.png'},
            {'size':(20,20)})

        plot_spectral(
            {'data':data2[col_idx].replace(np.nan, 0).as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data2',
                'filename':'dfH2.png'},
            {'size':(20,20)})


        # CALCULATE SD BY ROW
        # Subset by SD on data2. Cutoff = 1
        sd1 = data1[col_idx].std(axis=1)
        data1['sd'] = sd1.values
        plt.figure(figsize=(10,10))
        sns.distplot(sd1.values[~np.isnan(sd1.values)])
        plt.savefig('plots/dfH1_sdDistr.png')

        # Subset by SD on data2. Cutoff = 0.5
        sd2 = data2[col_idx].std(axis=1)
        data2['sd'] = sd2.values
        plt.figure(figsize=(10,10))
        sns.distplot(sd2.values[~np.isnan(sd2.values)])
        plt.savefig('plots/dfH2_sdDistr.png')

        print '\nShape after subsetting data1 on SD',data1[data1['sd']>1].shape
        print 'Shape after subsetting data2 on SD',data2[data2['sd']>0.5].shape


        # PLOT HEATMAPS SHOWING DATA AFTER SD-SUBSETTING
        # Subset by SD on data1. Cutoff = 1.0
        plot_spectral(
            {'data':data1[data1['sd'] > 0.5][col_idx].as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data1, SD > 0.5',
                'filename':'dfH1_sdSubset.png'},
            {'size':(20,20)})

        plot_spectral(
            {'data':data2[data2['sd'] > 0.5][col_idx].as_matrix(),
                'x':col_idx,'y':[]},
            {'title':'Data2, SD > 0.5',
                'filename':'dfH2_sdSubset.png'},
            {'size':(20,20)})

        # SAVE DATAFRAMES TO FILE
        data1[cols_to_save][data1['sd']>1].replace(np.nan, 0).to_csv(
            'data/dfH1_sdSubset.tsv', sep='\t', index=False)
        data2[cols_to_save][data2['sd']>0.5].replace(np.nan, 0).to_csv(
            'data/dfH2_sdSubset.tsv', sep='\t', index=False)

    if HL_analysis:
        dfHL = pd.merge(dfL, dfH, on='sym', how='outer')
        # Only keep data if represented across all variables.
        # Discards ~1000 genes
        dfHL = dfHL.dropna()

        prefix = 'ln_'

        row_idx = dfHL['sym']
        col_idx = [col for col in dfHL.columns if 'ln_' in col]

        data = dfHL[col_idx].as_matrix()
        data[data==0] = 0.001 # Clustering does not accept zeros

        # Bi- and Co-clustering
        n_clusters = len(col_idx)
        n_clusters2 = 97
        bc_model, bc_fit, bc_idx = Spectral_BiClustering(
            data, row_idx, col_idx, n_clusters)
        cc_model, cc_fit, cc_idx = Spectral_CoClustering(
            data, row_idx, col_idx, n_clusters2)

        # HEATMAPS:
        # Heatmap of original data
        plot_spectral(
            {'data':data,'x':col_idx,'y':row_idx},
            {'title':'HL: Original Data',
                'filename':prefix+'HL_original_data.png'},
            {'size':(20,20)})

        # Heatmap of biclustered data
        plot_spectral(
            {'data':bc_fit,'x':bc_idx[1],'y':bc_idx[0]},
            {'title':"HL: After biclustering; rearranged to show biclusters",
                'filename':prefix+'HL_bipartite_clustering.png'},
            {'size':(20,20)})

        # Heatmap of coclustered data
        plot_spectral(
            {'data':cc_fit,'x':cc_idx[1],'y':cc_idx[0]},
            {'title':"HL: After coclustering; rearranged to show coclusters",
                'filename':prefix+'HL_co_clustering.png'},
            {'size':(20,20)})

        # Compare bi-/co-clustered networks by Jaccard index
        jaccard_consensus(bc_model.biclusters_, cc_model.biclusters_)

    if AL_analysis:
        # Log-transformed columns
        prefix = 'ln_'

        row_idx = dfL['sym']
        col_idx = [col for col in dfL.columns if 'ln_' in col]

        data = dfL[col_idx].as_matrix()
        data[data==0] = 0.001 # Clustering does not accept zeros

        # Comparison of control experiments
        xy_linregress_plot(
            dfL[['C12.rna', 'C3.rna']],
            'C12.rna',
            'C3.rna',
            'Loyaza C12 vs. C3 RNA-Seq',
            'loayza_C12,3_rna.png')
        xy_linregress_plot(
            dfL[['C12.rp', 'C3.rp']],
            'C12.rp',
            'C3.rp',
            'Loyaza C12 vs. C3 Ribo-Seq',
            'loayza_C12,3_rp.png')

        # Bi- and Co-clustering
        bc_model, bc_fit, bc_idx = Spectral_BiClustering(
            data, row_idx, col_idx, len(col_idx))
        cc_model, cc_fit, cc_idx = Spectral_CoClustering(
            data, row_idx, col_idx, len(col_idx))

        # HEATMAPS:
        # Heatmap of original data
        plot_spectral(
            {'data':data,'x':col_idx,'y':row_idx},
            {'title':'L: Original Data',
                'filename':prefix+'L_original_data.png'},
            {'size':(10,20)})

        # Heatmap of biclustered data
        plot_spectral(
            {'data':bc_fit,'x':bc_idx[1],'y':bc_idx[0]},
            {'title':"L: After biclustering; rearranged to show biclusters",
                'filename':prefix+'L_bipartite_clustering.png'},
            {'size':(10,20)})

        # Heatmap of coclustered data
        plot_spectral(
            {'data':cc_fit,'x':cc_idx[1],'y':cc_idx[0]},
            {'title':"L: After coclustering; rearranged to show coclusters",
                'filename':prefix+'L_co_clustering.png'},
            {'size':(10,20)})

        # Compare bi-/co-clustered networks by Jaccard index
        jaccard_consensus(bc_model.biclusters_, cc_model.biclusters_)


