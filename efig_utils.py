#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import csv
import os
import sys
import re
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

import numpy as np
import scipy
from scipy.stats import mannwhitneyu as u
from scipy.stats import ttest_ind as t
from scipy import stats
from scipy.stats import norm
import math
from sklearn.metrics import r2_score

from itertools import cycle
from sklearn import svm, datasets
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from sklearn.metrics import roc_auc_score
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import PoissonRegressor
from sklearn.decomposition import PCA
from sklearn.decomposition import IncrementalPCA
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
import pickle
import statsmodels.api as sm

def sem(x,n):
    if n == 0: return 0.0
    else:
        p=float(x)/n
        return np.sqrt(p*abs((1-p))/float(n))

def ciOfOdds(x1,n1,x2,n2):
    import numpy as np
    import math
    a = x1
    b = n1-x1
    c = x2
    d = n2-x2
    if 0 in [a,b,c,d]: 
        a +=0.5
        b +=0.5
        c +=0.5
        d +=0.5
    odds = (a/b)/(c/d)
    log_odds = np.log(odds)
    se_log_odds = math.sqrt(1/a+1/b+1/c+1/d)
    upper_log_odds = log_odds + 1.96*se_log_odds
    lower_log_odds = log_odds - 1.96*se_log_odds
    return [odds,math.exp(lower_log_odds),math.exp(upper_log_odds)]


def download_fig_table(file_name):
    if not os.path.isdir('fig_tables'): os.system('mkdir fig_tables')
    if not os.path.isfile('fig_tables/{0}'.format(file_name)): 
        os.system('gsutil cp gs://gnomad-nc-constraint-v31-paper/fig_tables/{0} fig_tables/'.format(file_name))


# efig. 1a
def plt_po_mu(savefig):

    download_fig_table('fitted_po_by_context_methyl.txt')
    df_mr = pd.read_csv('fig_tables/fitted_po_by_context_methyl.txt',sep='\t',
                       index_col = ['context','ref','alt','methylation_level'])

    po = dict(zip(df_mr.index, df_mr.proportion_observed))
    mu = dict(zip(df_mr.index, df_mr.mu))
    po_type = dict(zip(df_mr.index, df_mr.variant_type))
    colors = sns.color_palette('Blues',40)[10:]
    d_color = {}
    for mut in po_type:
        if po_type[mut] == 'CpG':
            d_color[mut] = colors[int(mut[-1])]
        elif po_type[mut] == 'transversion': d_color[mut] = '#fdae61'
        elif po_type[mut] == 'non-CpG transition': d_color[mut] = '#66c2a5'

    plt.clf()
    plt.figure(figsize=(6,4))

    ks = po.keys()
    xp = [mu[k] for k in ks]
    yp = [po[k] for k in ks]
    c = [d_color[k] for k in ks]

    plt.scatter(xp, yp, color = c, alpha=0.8)
    plt.plot(sorted(df_mr['mu']), sorted(df_mr['fitted_po']), color='#969696', linestyle='dashed')

    # add legend
    k = ('CCG', 'C', 'T', 10)
    p1 = plt.scatter([mu[k]], [po[k]], color = d_color[k], alpha=0.8)

    vt = [i for i in po_type if po_type[i]=='non-CpG transition']
    l = [(i,po[i]) for i in vt]
    l.sort(key=lambda x:x[1])
    k = l[-1][0]
    p2 = plt.scatter([mu[k]], [po[k]], color = d_color[k], alpha=0.8)

    vt = [i for i in po_type if po_type[i]=='transversion']
    l = [(i,po[i]) for i in vt]
    l.sort(key=lambda x:x[1])
    k = l[0][0]
    p3 = plt.scatter([mu[k]], [po[k]], color = d_color[k], alpha=0.8)
    plt.legend([p1,p2,p3],['CpG transition','Non-CpG transition','Transversion'],loc='lower right')

    plt.xlabel('mu',fontsize=12.)
    plt.ylabel('Proportion observed',fontsize=12.)

    plt.ticklabel_format(axis = 'x', scilimits=[-4,6])

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 1b
def plt_po_mu_chrx(savefig):

    download_fig_table('fitted_po_by_context_methyl_chrx.txt')
    df_mr = pd.read_csv('fig_tables/fitted_po_by_context_methyl_chrx.txt',sep='\t',
                       index_col = ['context','ref','alt','methylation_level'])

    po = dict(zip(df_mr.index, df_mr.proportion_observed))
    mu = dict(zip(df_mr.index, df_mr.mu))
    po_type = dict(zip(df_mr.index, df_mr.variant_type))
    colors = sns.color_palette('Blues',40)[10:]
    d_color = {}
    for mut in po_type:
        if po_type[mut] == 'CpG':
            d_color[mut] = colors[int(mut[-1])]
        elif po_type[mut] == 'transversion': d_color[mut] = '#fdae61'
        elif po_type[mut] == 'non-CpG transition': d_color[mut] = '#66c2a5'

    plt.clf()
    plt.figure(figsize=(6,4))

    ks = po.keys()
    xp = [mu[k] for k in ks]
    yp = [po[k] for k in ks]
    c = [d_color[k] for k in ks]

    plt.scatter(xp, yp, color = c, alpha=0.8)
    plt.plot(sorted(df_mr['mu']), sorted(df_mr['fitted_po']), color='#969696', linestyle='dashed')
    # add legend
    k = ('CCG', 'C', 'T', 10)
    p1 = plt.scatter([mu[k]], [po[k]], color = d_color[k], alpha=0.8)

    vt = [i for i in po_type if po_type[i]=='non-CpG transition']
    l = [(i,po[i]) for i in vt]
    l.sort(key=lambda x:x[1])
    k = l[-1][0]
    p2 = plt.scatter([mu[k]], [po[k]], color = d_color[k], alpha=0.8)

    vt = [i for i in po_type if po_type[i]=='transversion']
    l = [(i,po[i]) for i in vt]
    l.sort(key=lambda x:x[1])
    k = l[0][0]
    p3 = plt.scatter([mu[k]], [po[k]], color = d_color[k], alpha=0.8)
    plt.legend([p1,p2,p3],['CpG transition','Non-CpG transition','Transversion'],loc='lower right')

    plt.xlabel('mu',fontsize=12.)
    plt.ylabel('Proportion observed',fontsize=12.)
    plt.ticklabel_format(axis = 'x', scilimits=[-4,6])

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 1c
def plt_genomic_ft_sel(savefig):
    
    download_fig_table('genomic_features13_sel.annot.txt')
    download_fig_table('mutation_rate_by_context_methyl.txt')
    dft = pd.read_csv('fig_tables/genomic_features13_sel.annot.txt',sep='\t')
    dft = dft.replace(np.nan,'')
    d_coef = dft.groupby(['feature','window'])[['context','coef']].apply(
        lambda x: x.set_index('context').to_dict(orient='index')).to_dict()

    d_coef_l = {}
    contexts = set([line.strip().split('\t')[0] for line in 
                     open('fig_tables/mutation_rate_by_context_methyl.txt').readlines()[1:]])
    contexts = sorted(list(contexts))
    ft_cols = ['GC_content', 'LCR', 'SINE', 'LINE', 'dist2telo', 'dist2cent',  
            'recomb_male', 'recomb_female', 
            'met_sperm','CpG_island','Nucleosome','cDNM_maternal_05M', 'cDNM_paternal_05M',]
    d_rename = {
     'cDNM_paternal_05M':'cDNM\npaternal',
     'GC_content':'GC content',
     'LINE':'LINE',
     'recomb_female':'Recomb rate\nfemale',
     'LCR':'LCR',
     'dist2cent':'Dist to\ncentromere',
     'CpG_island':'CpG island',
     'SINE':'SINE',
     'met_sperm':'Methyl',
     'dist2telo':'Dist to\ntelomere',
     'recomb_male':'Recomb rate\nmale',
     'Nucleosome':'Nucleosome',
     'cDNM_maternal_05M': 'cDNM\nmaternal',
    }
    
    ws = ['1k','10k','100k','1M']
    
    for ft in ft_cols:
        coef_l = []
        for window in ws:
            coef_l.append([d_coef[(ft,window)][i]['coef'] for i in contexts])
        d_coef_l[ft] = coef_l
    d_coef_label = dft.groupby(['feature','window'])[['context','label']].apply(
        lambda x: x.set_index('context').to_dict(orient='index')).to_dict()
    d_coef_l_label = {}
    for ft in ft_cols:
        coef_l_label = []
        for window in ws:
            coef_l_label.append([d_coef_label[(ft,window)][i]['label'] for i in contexts])
        d_coef_l_label[ft] = coef_l_label

    plt.clf()
    fig, axes = plt.subplots(nrows=1, ncols=len(ft_cols),  figsize=(15,len(contexts)*0.4),
                             sharex='col',sharey='row')
    plt.subplots_adjust(
        wspace=0.1, 
        hspace=0.1,)
    for i in range(0,len(ft_cols)):
        sns.heatmap(
            np.transpose(np.array(d_coef_l[ft_cols[i]])),
            square = True,
            vmin = -0.1, 
            vmax = 0.1,
            center = 0,
            cmap=sns.diverging_palette(220, 20, as_cmap=True),
            cbar=False,
            fmt = 's',
            annot = np.transpose(np.array(d_coef_l_label[ft_cols[i]])),
            annot_kws = {'color':'#737373'},
            ax = axes[i]
                    )
        if i==0: 
            axes[i].tick_params(axis='both',top=False,right=False,bottom=True,labelsize=10,rotation=0)
        else: 
            axes[i].tick_params(axis='both',top=False,right=False,bottom=True,left=False)
        axes[i].set_xticks([x+0.5 for x in range(0,len(ws))], labels=ws, rotation=90,fontsize=9. )
        axes[i].set_xlabel(d_rename[ft_cols[i]], rotation=0,fontsize=10.)
    axes[0].set_yticks([y+0.5 for y in range(0,len(contexts))],contexts,rotation=0,fontsize=10.)
    print ()

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 1d
def plt_scatter_oe_gnocchi(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    
    plt.clf()
    plt.figure(figsize=(6, 6))

    x = df_z['expected']
    y = df_z['observed']
    c = df_z['z']

    plt.scatter(x=x, y=y, c=c, cmap=sns.diverging_palette(220, 20, as_cmap=True), alpha=0.5, vmin=-10, vmax=10)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gnocchi score')
    
    markz = 4
    df_mark = df_z[ (df_z['z']>=markz-0.01) & (df_z['z']<=markz+0.01)]
    x_ = df_mark['expected']
    y_ = df_mark['observed']
    plt.plot(sorted(x_), sorted(y_), label = 'Z={0} (99th percentile)'.format(markz), 
            color='#fb6a4a', linestyle='dashed')
    plt.legend(loc='lower right')

    plt.ylabel('Observed number of variants',fontsize = 12.,)
    plt.xlabel('Expected number of variants',fontsize = 12.,)

    plt.xticks(np.arange(50, max(x)+1, 50))
    plt.yticks(np.arange(50, max(y)+1, 50))

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        savefig = savefig.replace('.pdf','.png')
        plt.savefig(savefig, bbox_inches='tight', dpi=600) 


# efig. 1e
def plt_scatter_oe_gnocchi_chrx(savefig):
    
    download_fig_table('constraint_z_genome_1kb_chrx.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb_chrx.annot.txt',sep='\t')
    
    plt.clf()
    plt.figure(figsize=(6, 6))
    
    x = df_z['expected']
    y = df_z['observed']
    c = df_z['z']
    
    plt.scatter(x=x, y=y, c=c, cmap=sns.diverging_palette(220, 20, as_cmap=True), alpha=0.5, vmin=-10, vmax=10)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Gnocchi score')

    markz = 4.7
    df_mark = df_z[ (df_z['z']>=markz-0.01) & (df_z['z']<=markz+0.01)]
    x_ = df_mark['expected']
    y_ = df_mark['observed']
    plt.plot(sorted(x_), sorted(y_), label = 'Z={0} (99th percentile)'.format(markz), 
             color='#fb6a4a', linestyle='dashed')
    plt.legend(loc='lower right')

    plt.ylabel('Observed number of variants',fontsize = 12.,)
    plt.xlabel('Expected number of variants',fontsize = 12.,)

    plt.xticks(np.arange(50, max(x)+1, 50))
    plt.yticks(np.arange(50, max(y)+1, 50))

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        savefig = savefig.replace('.pdf','.png')
        plt.savefig(savefig, bbox_inches='tight', dpi=600) 


# efig. 2a
def plt_prop_gnocchi4_cd(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    cut_bins = np.arange(0,1.1,.1)
    df_z['coding_prop_pct'] = pd.cut(df_z['coding_prop'], bins=cut_bins, 
                               labels=range(1,11))
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')
    df_z_exonic = df_z_exonic[df_z_exonic['pass_qc']]
    
    cutoff = 4
    x0, n0 = len(df_z[ (df_z['coding_prop']==0) & (df_z['z']>=cutoff) ]), len(df_z[df_z['coding_prop']==0])
    x0, n0 = 0,0
    x0 += len(df_z[ (df_z['coding_prop_pct']==1) & (df_z['z']>=cutoff) ])
    n0 += len(df_z[df_z['coding_prop_pct']==1])

    counts = [(x0,n0)]
    fracs = [x0/n0]
    sems = [sem(x0,n0)]
    for i in range(2,11):
        x, n = len(df_z[ (df_z['coding_prop_pct']==i) & (df_z['z']>=cutoff) ]), len(df_z[df_z['coding_prop_pct']==i])
        if i==10:
            x2, n2 = len(df_z_exonic[ df_z_exonic['z']>=cutoff ]), len(df_z_exonic)
            x += x2
            n += n2           
        counts.append((x,n))
        fracs.append(x/n)
        sems.append(sem(x,n))

    plt.clf()
    plt.figure(figsize=(6,4))
    colors = sns.cubehelix_palette(len(fracs))
    pos = list(range(0,len(fracs)))

    plt.bar(x = pos, height = [i*100 for i in fracs], yerr = [i*100 for i in sems],
            align='center',width=0.6,
            color=sns.cubehelix_palette(8)[2],
            alpha=.7,)

    plt.ylabel('Proportion of windows with\nGnocchi{1}{0} (%)'.format(cutoff,u'$\geq$'), fontsize=12)
    plt.xlabel('Percentage of coding bases (%)', fontsize=12)
    ticks = range(0,len(fracs))
    labels = list(range(0,100,10))+['Exonic\nonlly']
    ticks = range(0,len(labels))
    ticks = [i-0.5 for i in ticks]
    plt.xticks(ticks=ticks, labels=labels,)
    
    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 2b
def plt_hist_freq_gnocchi_exonic(savefig):
    
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')

    z1 = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]['z'].to_list()
    z2 = df_z_exonic[df_z_exonic['pass_qc']]['z'].to_list()
    
    plt.clf()
    fig, ax1 = plt.subplots(figsize=(6,4))
    color1, color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], sns.cubehelix_palette(8)[-3]
    color1, color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], sns.cubehelix_palette(8)[2]
    label1, label2 = 'Non-coding','Exonic only'
    hist_bins = np.arange(-10,11,0.5)

    ax1.hist(z1, bins=hist_bins, edgecolor=None, color = color1, alpha=0.7,density=False, label = label1)
    ax1.hist([0], edgecolor=None, color = color2, alpha=0.7,density=False, label = label2)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.hist(z2, bins=hist_bins, edgecolor=None, color = color2,alpha=0.7,density=False, label= label2)

    ax1.legend(loc='upper left', fontsize = 12.)
    ax1.set_xlabel('Gnocchi',fontsize = 12.)   
    ax1.set_ylabel('Frequency',fontsize = 12.,  color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax2.set_ylabel('Frequency',fontsize = 12.,  color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    
    plt.axvline(x = np.median(z1),color=color1,linestyle='dashed')
    plt.axvline(x = np.median(z2),color=color2,linestyle='dashed')
    plt.xlim(-10,10)
    plt.xticks(range(-10,11),range(-10,11))
    sns.despine(left=False, right=False, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=True)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 2c
def plt_prop_gnocchi4_cd_titration(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')

    z_nc = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]['z'].to_list()
    z_exonic = df_z_exonic[df_z_exonic['pass_qc']]['z'].to_list()

    cutoff = 4
    x0,n0 = len([i for i in z_nc if i>=cutoff]), len(z_nc)
    counts = [(x0,n0)]
    fracs = [x0/n0]
    sems = [sem(x0,n0)]

    import random
    random.seed(42)
    for p in np.arange(0.1,1.1,.1):
        z1 = random.sample(z_exonic, k=int(len(z_exonic)*p))
        z2 = z_nc+z1
        x, n = len([i for i in z2 if i>=cutoff]), len(z2)
        counts.append((x,n))
        fracs.append(x/n)
        sems.append(sem(x,n))

    plt.clf()
    plt.figure(figsize=(6,4))

    color = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    color = sns.cubehelix_palette(8)[2]
    plt.errorbar(y = [i*100 for i in fracs], x = range(0,len(fracs)), yerr=[i*100 for i in sems], 
                 color=color, ecolor=color, ls='', marker='o',elinewidth=2,alpha=0.7,)
    plt.ylabel('Proportion of windows with\nGnocchi{1}{0} (%)'.format(cutoff,u'$\geq$'), fontsize=12)
    plt.xlabel('Random draws of exonic windows (%)\ninto the non-coding genome', fontsize=12)
    plt.xticks(range(0,len(fracs)), np.arange(0,110,10))
     
    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 2d
def plt_gnochhi_ptcl_nc_cd(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')

    z_nc = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]['z'].to_list()
    z_exonic = df_z_exonic[df_z_exonic['pass_qc']]['z'].to_list()
    fracs = []
    for pctl in np.arange(10,110,10):
        z_pctl = np.percentile(z_exonic,pctl)
        pctl_nc = stats.percentileofscore(z_nc,z_pctl)
        fracs.append(pctl_nc)

    plt.clf()
    plt.figure(figsize=(6,4))
    
    color = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    color = sns.cubehelix_palette(8)[2]
    color_ = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    plt.scatter(y = fracs, x = range(0,len(fracs)), color=color, alpha=0.7)
    
    plt.ylabel('Gnocchi percentile of\nnon-coding region', fontsize=12)
    plt.xlabel('Gnocchi percentile of coding exon', fontsize=12)
    plt.xticks(range(0,len(fracs)), np.arange(10,110,10))
    color=sns.cubehelix_palette(8)[-3]
    color=sns.cubehelix_palette(8)[2]
    plt.axvline(x = 8,color=color,linestyle='dashed',ymax=fracs[8]/100)
    plt.axhline(y = fracs[8],color=color_,linestyle='dashed',xmax=.9)
    plt.axvline(x = 4,color=color,linestyle='dashed',ymax=fracs[4]/100)
    plt.axhline(y = fracs[4],color=color_,linestyle='dashed',xmax=.5)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 3a
def plt_gnocchi_ptcl_annot_cd(savefig):
    
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    download_fig_table('constraint_z_*_stringed_1kb.txt')
   
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')
    z_exonic = df_z_exonic[df_z_exonic['pass_qc']]['z'].to_list()
 
    plt.clf()
    plt.figure(figsize=(6,4))
    cmap = sns.cubehelix_palette(start=.5, rot=-.5, )
    cmap2 = sns.cubehelix_palette()
    ann_color = {   
        'Promoter': cmap2[-2],
        'Enhancer':cmap2[-4],
        'miRNA':cmap[-4],
        'lncRNA': cmap[-5],
    }
    ann_label = {'cCRE-PLS':'Promoter','cCRE-dELS':'Enhancer',
                 'lncRNA':'lncRNA',
                 'miRNA':'miRNA',                
                }
    for annot in ['cCRE-PLS','cCRE-dELS','miRNA','lncRNA']:
        df_z = pd.read_csv('fig_tables/constraint_z_{0}_stringed_1kb.txt'.format(annot),sep='\t',index_col='element_id')
        z_annot = df_z[(df_z['pass_qc'])]['z'].to_list()
        fracs = []
        for pctl in np.arange(10,110,10):
            z_pctl = np.percentile(z_exonic,pctl)
            pctl_nc = stats.percentileofscore(z_annot,z_pctl)
            fracs.append(pctl_nc)
        fracs2 = []
        step=2
        for pctl in np.arange(10,100+step,step):
            z_pctl = np.percentile(z_exonic,pctl)
            pctl_nc = stats.percentileofscore(z_annot,z_pctl)
            fracs2.append(pctl_nc)    
        
        color = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
        color = sns.cubehelix_palette(8)[2]
        color_ = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
        ann = ann_label[annot]
        plt.scatter(y = fracs, x = [i/10 for i in np.arange(10,110,10)], color=ann_color[ann], alpha=0.7,label=ann)
        plt.plot([i/10 for i in np.arange(10,100+step,step)], fracs2, color=ann_color[ann], linestyle='dashed')

    plt.legend(fontsize=12)
    plt.ylabel('Gnocchi percentile of\nregulatory sequence'.format(annot), fontsize=12)
    plt.xlabel('Gnocchi percentile of coding exon', fontsize=12)
    plt.xticks(range(1,len(fracs)+1), np.arange(10,110,10))
    color=sns.cubehelix_palette(8)[-3]
    color=sns.cubehelix_palette(8)[2]
    plt.axvline(x = 5,color='#969696',linestyle='dashed',ymax=fracs[5]/100)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
        
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 3b
def plt_gnocchi_ptcl_annot_nc(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_*_stringed_1kb.txt')
   
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    z_nc = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]['z'].to_list()  
 
    plt.clf()
    plt.figure(figsize=(6,4))
    cmap = sns.cubehelix_palette(start=.5, rot=-.5, )
    cmap2 = sns.cubehelix_palette()
    ann_color = {   
        'Promoter': cmap2[-2],
        'Enhancer':cmap2[-4],
        'miRNA':cmap[-4],
        'lncRNA': cmap[-5],
    }
    ann_label = {'cCRE-PLS':'Promoter','cCRE-dELS':'Enhancer',
                 'lncRNA':'lncRNA',
                 'miRNA':'miRNA',                
                }
    for annot in ['cCRE-PLS','cCRE-dELS','miRNA','lncRNA']:
        df_z = pd.read_csv('fig_tables/constraint_z_{0}_stringed_1kb.txt'.format(annot),sep='\t',index_col='element_id')
        z_annot = df_z[(df_z['pass_qc'])]['z'].to_list()
        fracs = []
        for pctl in np.arange(10,110,10):
            z_pctl = np.percentile(z_nc,pctl)
            pctl_nc = stats.percentileofscore(z_annot,z_pctl)
            fracs.append(pctl_nc)
        fracs2 = []
        step=2
        for pctl in np.arange(10,100+step,step):
            z_pctl = np.percentile(z_nc,pctl)
            pctl_nc = stats.percentileofscore(z_annot,z_pctl)
            fracs2.append(pctl_nc)    
        
        color = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
        color = sns.cubehelix_palette(8)[2]
        color_ = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
        ann = ann_label[annot]
        plt.scatter(y = fracs, x = [i/10 for i in np.arange(10,110,10)], color=ann_color[ann], alpha=0.7,label=ann)
        plt.plot([i/10 for i in np.arange(10,100+step,step)], fracs2, color=ann_color[ann], linestyle='dashed')

    plt.legend(fontsize=12)
    plt.ylabel('Gnocchi percentile of\nregulatory sequence'.format(annot), fontsize=12)
    plt.xlabel('Gnocchi percentile of coding exon', fontsize=12)
    plt.xticks(range(1,len(fracs)+1), np.arange(10,110,10))
    color=sns.cubehelix_palette(8)[-3]
    color=sns.cubehelix_palette(8)[2]
    plt.axvline(x = 5,color='#969696',linestyle='dashed',ymax=fracs[5]/100)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
        
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 4a
def plt_enrichment_gwas_vs_ccre(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')    
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    dft = df_z[ (df_z['coding_prop']==0) & (df_z['pass_qc'])]
    dft['ENCODE cCRE'] = dft[['ENCODE cCRE-{0}'.format(ann) for ann in ['PLS','pELS','dELS']]].sum(axis=1)
    dft['GWAS'] = dft['GWAS Catalog repl (ext)'] | dft['GWAS fine-mapping']
    
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])   
    label = int((len(cut_bins)-1)/2)
    dft['z_bin'] = pd.cut(dft['z'], bins=cut_bins, labels=range((-1)*label,label))
    
    plt.clf()
    plt.figure(figsize=(6,4))
    ann1,ann2 = 'GWAS','ENCODE cCRE'
    ann_color = { 
        'GWAS':  sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
    }
    
    for exl in [False, True]:
        if exl: 
            dfp=dft[dft[ann2]==0]
        else: 
            dfp=dft
            dfp=dft[dft[ann2]>0]
        odds2plot = []
        ci2plot = []
        l1 = dfp[ann1].dropna().to_list() # baseline
        x1,n1 = len([i for i in l1 if i]),len(l1)

        bins = sorted(list(set(dfp['z_bin'])))
        for z in bins:
            l2 = dfp[dfp['z_bin']==z][ann1].dropna().to_list()
            x2,n2 = len([i for i in l2 if i]),len(l2)
            ci = ciOfOdds(x2,n2,x1,n1)
            odds2plot.append(ci[0])
            ci2plot.append([ci[0]-ci[1],ci[2]-ci[0]])

        rightshift = 0.1 if exl else -0.1
        alpha = 0.4 if exl else .7
        label = 'Outside cCREs' if exl else 'Within cCREs'
        fmt = '-' if exl else '--'
        color = ann_color[ann1]
        plt.errorbar(y = odds2plot, 
                     x = [i+0.5+rightshift for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o', ecolor=color, color=color,fmt=fmt,linestyle='', elinewidth=2, alpha=alpha,
                     label = label,
                    )

    plt.axhline(y = 1.0,color='#969696',linestyle='dashed')
    plt.rcParams['legend.title_fontsize'] = 12.
    plt.legend(loc='upper left', fontsize = 12., title='GWAS variants')
    plt.xticks(bins[1:],bins[1:])
    plt.xlabel('Gnocchi', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 4b
def plt_delta_pip(savefig):
    
    download_fig_table('UKBB.nc_constraint.updated_pip.nc_in_cs.tsv')  
    df_pip = pd.read_csv('fig_tables/UKBB.nc_constraint.updated_pip.nc_in_cs.tsv',sep='\t')
    df_pip['delta_pip'] = df_pip['susie.pip_updated']-df_pip['susie.pip']

    pip_=.8
    d_pip = 0.01
    dfp = df_pip[(df_pip['target']>=0)
           & (df_pip['susie.pip_updated']>=pip_) 
           & (df_pip['susie.pip']<pip_) 
           & (df_pip['delta_pip']>d_pip)
                ]
                
    plt.clf()
    plt.figure(figsize=(6,4))

    y = dfp['delta_pip']
    x = dfp['susie.pip_updated']
    plt.scatter(x=x,y=y,color=sns.cubehelix_palette(start=.5, rot=-.5, )[-2], alpha=0.4)

    plt.xlabel('PIP$_G$$_n$$_o$$_c$$_c$$_h$$_i$', fontsize=12)
    plt.ylabel('PIP$_G$$_n$$_o$$_c$$_c$$_h$$_i$ - PIP$_u$$_n$$_i$$_f$', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight')


# efig. 5a
def plt_comparison_roc(pos,neg,dist2exon,savefig):
    
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')

    scores = ['Orion','CDTS','gwRVIS','DR',"phastCons","phyloP","GERP"]
    scores.append('z')
    score_color = {
        'z':sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
        'Orion':'#33a02c','CDTS':'#fb9a99','DR':'#993404','gwRVIS':'#6baed6',   
        'phastCons':'#969696','phyloP':'#737373','GERP':'#bdbdbd'
    }
    score_label = {'z': 'Gnocchi'} 

    df_pos = {'gwas_catalog': 'comparisons_gwas_catalog_repl', 
              'gwas_fine-mapping':'comparisons_gwas_fine-mapping_pip09', 
              'gwas_fine-mapping_hc':'comparisons_gwas_fine-mapping_pip09_hc', 
              'clinvar_plp_hgmd':'comparisons_41588_2018_62_MOESM4_ESM_clinvar_plp',     
          }
    title_label = {'gwas_catalog':'GWAS Catalog','gwas_fine-mapping':'GWAS fine-mapping',
                   'gwas_fine-mapping_hc':'GWAS fine-mapping\n(high confidence)',
                   'clinvar_plp_hgmd':'Likely pathogenic'}
    df_neg = {
        'topmed_maf5':'comparisons_topmed_maf5.sampled.cov', 
        'topmed_mac1':'comparisons_topmed_mac1.sampled.cov', 
          }
    sampling = 10
    
    plt.clf()
    fig,ax = plt.subplots(1, figsize=(4.5,4.5))

    df_1 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_pos[pos]),sep='\t')
    df_1 = df_1[df_1['dist2exon']>=dist2exon]
    df_1['group'] = 1
    df_0 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_neg[neg]),sep='\t')
    df_0 = df_0[df_0['dist2exon']>=dist2exon]
    df_0['group'] = 0 
    if sampling: df_0 = df_0.sample(n=sampling*len(df_1), random_state=714)
    df_01 = pd.concat([df_1,df_0]).drop_duplicates(subset=['locus'])

    d_auc = {}
    d_auc_ = {}
    idx = 0
    for score in scores: 
        dfp_01_p = df_01[[score,'group']].dropna()
        y_true = dfp_01_p['group']
        y_probas = dfp_01_p[score]
        fpr, tpr, _ = roc_curve(y_true,  y_probas)
        auc = roc_auc_score(y_true, y_probas)
        label = score_label[score] if score in score_label else score
        plt.plot(fpr,tpr,color=score_color[score], alpha=0.7, 
                 label='{0} ({1})'.format(label,round(auc,3)))
        d_auc[idx] = auc
        d_auc_[score] = auc
        idx +=1 
    plt.plot([0, 1], [0, 1], linestyle='--',color='#bdbdbd')

    sorted_auc = [i[0] for i in sorted(d_auc.items(), key=lambda item: item[1])][::-1]
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[i] for i in sorted_auc]
    labels = [labels[i] for i in sorted_auc]

    plt.legend(handles,labels,bbox_to_anchor=(1, 0.75))
    plt.title(title_label[pos])
    plt.ylabel('True positive rate', fontsize=12)
    plt.xlabel('False positive rate', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 5b
def plt_comparison_auc_af(pos,savefig):
    
    def count_n_maf(df):
        d = {}
        df = df.dropna(subset=['AF_topmed'])
        d['mac0'] = len(df[ (df['AC_topmed'] == 0)])
        d['mac1'] = len(df[ (df['AC_topmed'] == 1)])
        d['maf001'] = len(df[ (df['AC_topmed'] > 1) & (df['AF_topmed'] <= 0.01/100)])
        d['maf001_01'] = len(df[ (df['AF_topmed'] > 0.01/100) & (df['AF_topmed'] <= 0.1/100)])
        d['maf01_1'] = len(df[ (df['AF_topmed'] > 0.1/100) & (df['AF_topmed'] <= 1/100)])
        d['maf1_5'] = len(df[ (df['AF_topmed'] > 1/100) & (df['AF_topmed'] <= 5/100)])
        d['maf5'] = len(df[ (df['AF_topmed'] > 5/100)])
        d2 = dict([k,d[k]/len(df)*100] for k in d)
        return d,d2
    
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')

    scores = ['Orion','CDTS','gwRVIS','DR']
    scores.append('z')
    score_color = {
        'z':sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
        'Orion':'#33a02c','CDTS':'#fb9a99','DR':'#993404','gwRVIS':'#6baed6',   
        'phastCons':'#969696','phyloP':'#737373','GERP':'#bdbdbd'
    }
    score_label = {'z': 'Gnocchi'} 

    df_pos = {'gwas_catalog': 'comparisons_gwas_catalog_repl', 
              'gwas_fine-mapping':'comparisons_gwas_fine-mapping_pip09', 
              'gwas_fine-mapping_hc':'comparisons_gwas_fine-mapping_pip09_hc', 
              'clinvar_plp_hgmd':'comparisons_likely_pathogenic_clinvar_hgmd',     
          }
    title_label = {'gwas_catalog':'GWAS Catalog','gwas_fine-mapping':'GWAS fine-mapping',
                   'gwas_fine-mapping_hc':'GWAS fine-mapping\n(high confidence)',
                   'clinvar_plp_hgmd':'Likely pathogenic'}

    df_neg = {
        'topmed_maf5': 'comparisons_topmed_maf5.sampled.cov', 
        'topmed_mac1': 'comparisons_topmed_mac1.sampled.cov', 
        'topmed_maf001': 'comparisons_topmed_maf001.sampled.cov', 
        'topmed_maf001_01': 'comparisons_topmed_maf001_01.sampled.cov', 
        'topmed_maf01_1': 'comparisons_topmed_maf01_1.sampled.cov', 
        'topmed_maf1_5': 'comparisons_topmed_maf1_5.sampled.cov', 
          }

    sampling = 10

    df_1 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_pos[pos]),sep='\t')
    df_1['group'] = 1

    d_auc = {}
    for score in scores: 
        d_auc[score] = []
        for neg in ['topmed_mac1','topmed_maf001','topmed_maf001_01','topmed_maf01_1','topmed_maf1_5','topmed_maf5']:
            df_0 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_neg[neg]),sep='\t')
            df_0['group'] = 0
            if sampling: df_0 = df_0.sample(n=min(sampling*len(df_1),len(df_0)), random_state=714)

            df_01 = pd.concat([df_1,df_0]).drop_duplicates(subset=['locus'])

            dfp_01_p = df_01[[score,'group']].dropna()
            y_true = dfp_01_p['group']
            y_probas = dfp_01_p[score]
            fpr, tpr, _ = roc_curve(y_true,  y_probas)

            auc = roc_auc_score(y_true, y_probas)
            d_auc[score].append(auc)

    d_auc = {k: v for k, v in sorted(d_auc.items(), key=lambda item: item[1], reverse=True)}

    plt.clf()
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(4.5, 5), sharex = True,
                           gridspec_kw={'height_ratios': [1, 1.5]})
    for score in scores:
        if score not in d_auc: continue
        l = d_auc[score]
        label = score_label[score] if score in score_label else score
        axs[1].plot(range(0,len(l)),l, marker='o', linestyle='dashed',
                     label=label,color=score_color[score],alpha=0.7)

    count = list(count_n_maf(df_1)[1].values())[1:]
    axs[0].bar(x=range(0,len(count)),height=count,color='#bdbdbd',alpha=0.5)
    axs[0].set_ylabel('Proportion of positive\nvariant set (%)',fontsize = 12.)
    axs[0].set_ylim(0,100)

    axs[1].set_ylabel('AUC',fontsize = 12.)
    plt.xticks(range(0,len(count)),
               ['AC=1','(AC=1,\n0.01%]','(0.01%,\n0.1%]','(0.1%,\n1%]','(1%,\n5%]','>5%'])

    plt.xlabel('AF threshold for negative variant set',fontsize=12)
    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 6
def plt_comparison_roc_gnocchi(pos,neg,dist2exon,savefig):
    
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')

    scores = ['z','z_sliding100','z_trimer','z_heptamer',]
    score_label = {
        'z': 'Gnocchi',
        'z_sliding100':'Gnocchi, sliding',
        'z_trimer':'Trimer-only',
        'z_heptamer':'Heptamer-only',
                  }   
    c = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    c2 = sns.cubehelix_palette(start=.5, rot=-.5, )[-3]
    c3 = sns.cubehelix_palette(start=.5, rot=-.5, )[-4]
    c4 = sns.cubehelix_palette(start=.5, rot=-.5, )[-5]
    score_color = {'z':c,'z_sliding100':c2,'z_trimer':c3,'z_heptamer':c4,}

    df_pos = {'gwas_catalog': 'comparisons_gwas_catalog_repl', 
              'gwas_fine-mapping':'comparisons_gwas_fine-mapping_pip09', 
              'gwas_fine-mapping_hc':'comparisons_gwas_fine-mapping_pip09_hc', 
              'clinvar_plp_hgmd':'comparisons_likely_pathogenic_clinvar_hgmd',     
          }
    title_label = {'gwas_catalog':'GWAS Catalog','gwas_fine-mapping':'GWAS fine-mapping',
                   'gwas_fine-mapping_hc':'GWAS fine-mapping\n(high confidence)',
                   'clinvar_plp_hgmd':'Likely pathogenic'}

    df_neg = {
        'topmed_maf5': 'comparisons_topmed_maf5.sampled.cov', 
        'topmed_mac1': 'comparisons_topmed_mac1.sampled.cov', 
          }
    sampling = 10

    plt.clf()
    fig,ax = plt.subplots(1, figsize=(4.5,4.5))

    df_1 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_pos[pos]),sep='\t')
    df_1 = df_1[df_1['dist2exon']>=dist2exon]
    df_1['group'] = 1
    df_0 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_neg[neg]),sep='\t')
    df_0 = df_0[df_0['dist2exon']>=dist2exon]
    df_0['group'] = 0 
    if sampling: df_0 = df_0.sample(n=sampling*len(df_1), random_state=714)

    df_01 = pd.concat([df_1,df_0]).drop_duplicates(subset=['locus'])

    d_auc = {}
    d_auc_ = {}
    idx = 0
    for score in scores: 

        dfp_01_p = df_01[[score,'group']].dropna()

        y_true = dfp_01_p['group']
        y_probas = dfp_01_p[score]
        fpr, tpr, _ = roc_curve(y_true,  y_probas)

        auc = roc_auc_score(y_true, y_probas)
        label = score_label[score] if score in score_label else score
        plt.plot(fpr,tpr,color=score_color[score], alpha=0.7, 
                 label='{0} ({1})'.format(label,round(auc,3)))
        d_auc[idx] = auc
        d_auc_[score] = auc
        idx +=1 

    plt.plot([0, 1], [0, 1], linestyle='--',color='#bdbdbd')

    sorted_auc = [i[0] for i in sorted(d_auc.items(), key=lambda item: item[1])][::-1]
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[i] for i in sorted_auc]
    labels = [labels[i] for i in sorted_auc]

    ax.legend(handles,labels,loc="lower right")
    plt.title(title_label[pos])

    plt.ylabel('True positive rate', fontsize=12)
    plt.xlabel('False positive rate', fontsize=12)


    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 7
def plt_score_corr(scores,savefig):
    
    download_fig_table('genome_1kb.scores.txt')
    df_scores = pd.read_csv('fig_tables/genome_1kb.scores.txt',sep='\t')
    df_scores = df_scores.rename(columns={"Constraint Z":"Gnocchi"})

    plt.clf()
    scale = 0.7
    plt.figure(figsize=(6,6))

    corr = 'spearman'
    with sns.axes_style('white'):
        ax = sns.heatmap(
                         df_scores[scores].corr(method=corr),
                         square = True,
                         vmin = -1., vmax = 1.,
                         center = 0,
                         cmap=sns.diverging_palette(220, 20, as_cmap=True),
                         cbar=True,cbar_kws={'label': '{0} correlation'.format(corr.capitalize()),
                                            'shrink': 0.6},
                         linewidths = 0.2, linecolor= 'lightgrey',
                         annot = True,)

        plt.yticks(rotation=0,fontsize=12.)
        plt.xticks(rotation=90,fontsize=12.)
        plt.tight_layout()

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
    
# efig. 8a-b
def plt_power_depl(ws,savefig):
    
    download_fig_table('powered_z*_log-log.nc.txt')
    z = 4
    dfp_list = []
    for depl in list(np.arange(0.1,0.75,0.1)) + [0.2, 0.42, 0.48] + [0.26, 0.56, 0.64]:
        min_exp = math.ceil(z*z/depl/depl) 
        dft = pd.read_table('fig_tables/powered_z{0}-depl{1}-exp{2}_log-log.nc.txt'.format(
            z,'{:.2f}'.format(depl),min_exp))
        dft['depl'] = round(depl,3)
        dfp_list.append(dft)    
    dfp = pd.concat(dfp_list).drop_duplicates()

    # add in exon benchmarks
    from scipy.interpolate import make_interp_spline

    plt.clf()
    logscale=True

    dfp_w = dfp[dfp['Window_size'] == ws]
    dfp_w = dfp_w.sort_values(by='Sample_size')
    step = 0.1
    depl_min, depl_max = 0.1, 0.5
    l = np.arange(depl_min, depl_max+step,step)
    ### exon 800-1200bp
    if ws == 1000: 
        l_e = [0.2, 0.42, 0.48]
        label_sfx = {0.2: '(50th pctl exon)', 0.42: '(90th pctl exon)', 0.48: '(95th pctl exon)'}
    ### exon 80-120bp
    if ws == 100: 
        l_e = [0.26, 0.56, 0.64]
        label_sfx = {0.26: '(50th pctl exon)', 0.56: '(90th pctl exon)', 0.64: '(95th pctl exon)'}
    l = sorted(set(list(l)+l_e))
    l = [round(i,3) for i in l]

    colors = sns.cubehelix_palette(len(l), start=.5, rot=-.5,)

    for depl in l:

        x = dfp_w[dfp_w['depl']==depl]['Sample_size'].to_list()
        y = dfp_w[dfp_w['depl']==depl]['Percent'].to_list()

        i = l.index(depl)
        color = colors[i]
        lw = 3
        label = '{0}%'.format(round(depl*100))
        ls = 'solid'
        if depl in label_sfx:
            label += '\n'+label_sfx[depl]
            ls = 'dotted'
        plt.plot(x, y, color=color, lw=lw, label=label, ls=ls )

    plt.legend(
        title='Depletion of variation',
        fontsize = 11.,bbox_to_anchor=(1, 1), )

    plt.yticks([0,.2,.4,.6,.8, 1.0], ['0','20','40','60','80','100'])
    if logscale: plt.xscale('log')
    plt.xlabel('Sample size required',fontsize = 12.,)
    if ws == 1000: plt.ylabel('Percent of 1kb windows powered\nto detect constraint (Gnocchi{0}4)'.format(u'$\geq$'),fontsize = 12.,)
    if ws == 100: plt.ylabel('Percent of 100bp windows powered\nto detect constraint (Gnocchi{0}4)'.format(u'$\geq$'), fontsize = 12.,)

    n_v3 = 76156
    plt.axvline(x = n_v3, color='#969696', linestyle='dashed',lw=2.0)
    plt.axhline(y = .8, color='#fcbba1',)
    
    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 8c
def plt_comparison_auc_ws(savefig):
    
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')
    
    ws = ['100bp','500bp','2kb','3kb']
    scores = ['z']
    scores += ['z_{0}'.format(w) for w in ws]

    df_pos = {'gwas_catalog': 'comparisons_gwas_catalog_repl', 
              'gwas_fine-mapping':'comparisons_gwas_fine-mapping_pip09', 
              'gwas_fine-mapping_hc':'comparisons_gwas_fine-mapping_pip09_hc', 
              'clinvar_plp_hgmd':'comparisons_likely_pathogenic_clinvar_hgmd',     
          }
    title_label = {'gwas_catalog':'GWAS Catalog','gwas_fine-mapping':'GWAS fine-mapping',
                   'gwas_fine-mapping_hc':'GWAS fine-mapping\n(high confidence)',
                   'clinvar_plp_hgmd':'Likely pathogenic'}

    df_neg = {
        'topmed_maf5': 'comparisons_topmed_maf5.sampled.cov', 
        'topmed_mac1': 'comparisons_topmed_mac1.sampled.cov', 
          }

    sampling = 10

    neg = 'topmed_maf5'
    d_dat_auc = {}
    for pos in ['gwas_catalog','gwas_fine-mapping','gwas_fine-mapping_hc','clinvar_plp_hgmd']:
        if pos == 'clinvar_plp_hgmd': neg = 'topmed_mac1'
        d_dat_auc[pos] = {}
        df_1 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_pos[pos]),sep='\t')
        df_1['group'] = 1
        df_0 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_neg[neg]),sep='\t')
        df_0['group'] = 0 
        if sampling: df_0 = df_0.sample(n=sampling*len(df_1), random_state=714)

        df_01 = pd.concat([df_1,df_0]).drop_duplicates(subset=['locus'])

        for score in scores: 
            dfp_01_p = df_01[[score,'group']].dropna()
            y_true = dfp_01_p['group']
            y_probas = dfp_01_p[score]
            fpr, tpr, _ = roc_curve(y_true,  y_probas)
            auc = roc_auc_score(y_true, y_probas)
            d_dat_auc[pos][score] = auc
    d_ws = {'100':'z_100bp','500':'z_500bp','1000':'z','2000':'z_2kb','3000':'z_3kb',}
    d_color = {'clinvar_plp_hgmd':sns.cubehelix_palette(8)[2],
               'gwas_catalog':sns.cubehelix_palette(start=.5, rot=-.5, )[-4],
               'gwas_fine-mapping':sns.cubehelix_palette(start=.5, rot=-.5, )[-3],
               'gwas_fine-mapping_hc':sns.cubehelix_palette(start=.5, rot=-.5, )[-2]}
    d_label = title_label

    plt.clf()
    plt.figure(figsize=(4.5,4.5))
    for dat in d_dat_auc:
        y = []
        for ws in d_ws:
            score = d_ws[ws]
            y.append(d_dat_auc[dat][score])  

        plt.plot(range(0,len(y)),y, marker='o', linestyle='dashed',
                 label=d_label[dat],color=d_color[dat],alpha=0.7)

    plt.xticks(range(0,len(d_ws)),d_ws.keys())
    plt.legend(loc='lower right')
    plt.ylim(0.5,0.75)
    locs, labels = plt.xticks()
    plt.ylabel('AUC',fontsize=12)
    plt.xlabel('Window size (bp)',fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 8d
def plt_comparison_auc_pop(savefig):
    
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')
    pops = ['global','nfe']

    scores = ['z']
    scores += ['z_{0}'.format(p) for p in pops]

    df_pos = {'gwas_catalog': 'comparisons_gwas_catalog_repl', 
              'gwas_fine-mapping':'comparisons_gwas_fine-mapping_pip09', 
              'gwas_fine-mapping_hc':'comparisons_gwas_fine-mapping_pip09_hc', 
              'clinvar_plp_hgmd':'comparisons_likely_pathogenic_clinvar_hgmd',     
          }
    title_label = {'gwas_catalog':'GWAS Catalog','gwas_fine-mapping':'GWAS fine-mapping',
                   'gwas_fine-mapping_hc':'GWAS fine-mapping\n(high confidence)',
                   'clinvar_plp_hgmd':'Likely pathogenic'}

    df_neg = {
        'topmed_maf5': 'comparisons_topmed_maf5.sampled.cov', 
        'topmed_mac1': 'comparisons_topmed_mac1.sampled.cov', 
          }

    sampling = 10

    neg = 'topmed_maf5'
    d_dat_auc = {}
    for pos in ['gwas_catalog','gwas_fine-mapping','gwas_fine-mapping_hc','clinvar_plp_hgmd']:
        if pos == 'clinvar_plp_hgmd': neg = 'topmed_mac1'
        d_dat_auc[pos] = {}
        df_1 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_pos[pos]),sep='\t')
        df_1['group'] = 1
        df_0 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_neg[neg]),sep='\t')
        df_0['group'] = 0 
        if sampling: df_0 = df_0.sample(n=sampling*len(df_1), random_state=714)

        df_01 = pd.concat([df_1,df_0]).drop_duplicates(subset=['locus'])

        for score in scores: 
            dfp_01_p = df_01[[score,'group']].dropna()
            y_true = dfp_01_p['group']
            y_probas = dfp_01_p[score]
            fpr, tpr, _ = roc_curve(y_true,  y_probas)
            auc = roc_auc_score(y_true, y_probas)
            d_dat_auc[pos][score] = auc
    c = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    c2 = sns.cubehelix_palette(start=.5, rot=-.5, )[-3]
    c3 = sns.cubehelix_palette(start=.5, rot=-.5, )[-4]
    c4 = sns.cubehelix_palette(start=.5, rot=-.5, )[-5]
    d_color = {'z':c,
               'z_global':c3,
               'z_nfe':c2,}
    d_label = {'z':'gnomAD N=76,156',
               'z_global':'gnomAD downsampled N=34,029',
               'z_nfe':'gnomAD NFE N=34,029',}
    d_label_x = {'gwas_catalog':'GWAS\nCatalog','gwas_fine-mapping':'GWAS\nfine-mapping',
                 'gwas_fine-mapping_hc':'GWAS\nfine-mapping\n(high confidence)',
                 'clinvar_plp_hgmd':'Likely\npathogenic'}
    zs = ['z','z_nfe','z_global']
    
    plt.clf()
    plt.figure(figsize=(4.5,4.5))
    for z in zs:
        plt.scatter(x=range(0,len(d_dat_auc)),y=[d_dat_auc[dat][z] for dat in d_dat_auc],
                    color=d_color[z],alpha=0.7,label=d_label[z])

    plt.legend(loc='upper left')
    plt.ylim(0.6,0.8)
    plt.ylabel('AUC',fontsize=12)
    plt.xticks(range(0,len(d_dat_auc)),[d_label_x[i] for i in d_dat_auc],rotation=0)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 




