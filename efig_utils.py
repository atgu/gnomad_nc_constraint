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
def plt_pred_dnms_ft(savefig):    
    
    download_fig_table('rf_f18_predicted_dnms_1M.txt')
    download_fig_table('rf_f18_feature_importance.txt')

    df_dmn = pd.read_csv('fig_tables/rf_f18_predicted_dnms_1M.txt',sep='\t')
    df_imp = pd.read_csv('fig_tables/rf_f18_feature_importance.txt',sep='\t')
    
    plt.clf()
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(4, 4+len(df_imp)*0.4),
                           gridspec_kw={'height_ratios': [4, len(df_imp)*0.4],})
    # scatter plot
    y = df_dmn['observed']
    x = df_dmn['predicted']
    axs[0].scatter(x, y, color = '#74add1', alpha=0.3)
    m1 = min(min(x),min(y))
    m2 = max(max(x),max(y))
    axs[0].set_ylabel('Observed DNMs',fontsize=12.)
    axs[0].set_xlabel('Predicted DNMs',fontsize=12.)
    axs[0].plot([m1, m2], [m1, m2], color='grey', linestyle='--')
    # bar plot
    axs[1].barh(range(len(df_imp)), df_imp['importance'][::-1],
            color='#74add1', align='center')
    axs[1].set_xlabel('Feature importance',fontsize = 12.)
    axs[1].set_yticks(range(len(df_imp)))
    axs[1].set_yticklabels(list(df_imp['feature'][::-1]))
    axs[1].set_ylim(-0.5,17.5)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 

# efig. 1d
def plt_scatter_oe_z(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    
    plt.clf()
    plt.figure(figsize=(6, 6))

    x = df_z['expected']
    y = df_z['observed']
    c = df_z['z']

    plt.scatter(x=x, y=y, c=c, cmap=sns.diverging_palette(220, 20, as_cmap=True), alpha=0.5, vmin=-10, vmax=10)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Constraint Z score')
    
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
def plt_scatter_oe_z_chrx(savefig):
    
    download_fig_table('constraint_z_genome_1kb_chrx.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb_chrx.annot.txt',sep='\t')
    
    plt.clf()
    plt.figure(figsize=(6, 6))
    
    x = df_z['expected']
    y = df_z['observed']
    c = df_z['z']
    
    plt.scatter(x=x, y=y, c=c, cmap=sns.diverging_palette(220, 20, as_cmap=True), alpha=0.5, vmin=-10, vmax=10)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Constraint Z score')

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
def plt_hist_freq_z_exonic(savefig):
    
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')
    z1 = df_z[df_z['coding_prop']==0]['z'].to_list()
    z2 = df_z_exonic['z'].to_list()
    
    plt.clf()
    fig, ax1 = plt.subplots(figsize=(6,4))

    color1, color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], sns.cubehelix_palette(8)[-3]
    label1, label2 = 'Non-coding','Exonic only'
    hist_bins = np.arange(-10,11,1)

    ax1.hist(z1, bins=hist_bins, edgecolor=None, color = color1, alpha=0.7,density=False, label = label1)
    ax1.hist([0], edgecolor=None, color = color2, alpha=0.7,density=False, label = label2)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.hist(z2, bins=hist_bins, edgecolor=None, color = color2,alpha=0.7,density=False, label= label2)

    ax1.legend(loc='upper left', fontsize = 12.)
    ax1.set_xlabel('Constraint Z',fontsize = 12.)   
    ax1.set_ylabel('Frequency',fontsize = 12.,  color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax2.set_ylabel('Frequency',fontsize = 12.,  color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    
    plt.axvline(x = 0.0,color='#525252',linestyle='dashed')
    plt.xlim(-10,10)
    plt.xticks(range(-10,11),range(-10,11))
    sns.despine(left=False, right=False, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=True)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 

# efig. 2b
def plt_prop_z4_cd(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',
                       index_col='element_id')[['z','coding_prop']]
    cut_bins = np.arange(0,1.1,.1)
    df_z['coding_prop_pct'] = pd.cut(df_z['coding_prop'], bins=cut_bins, 
                               labels=range(1,11))
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')

    cutoff = 4
    x0, n0 = len(df_z[ (df_z['coding_prop']==0) & (df_z['z']>=cutoff) ]), len(df_z[df_z['coding_prop']==0])
    x0 += len(df_z[ (df_z['coding_prop_pct']==1) & (df_z['z']>=cutoff) ])
    n0 += len(df_z[df_z['coding_prop_pct']==1])

    counts = [(x0,n0)]
    fracs = [x0/n0]
    sems = [sem(x0,n0)]
    for i in range(2,11):
        x, n = len(df_z[ (df_z['coding_prop_pct']==i) & (df_z['z']>=cutoff) ]), len(df_z[df_z['coding_prop_pct']==i])
        counts.append((x,n))
        fracs.append(x/n)
        sems.append(sem(x,n))

    x, n = len(df_z_exonic[ df_z_exonic['z']>=cutoff ]), len(df_z_exonic)
    counts.append((x,n))
    fracs.append(x/n)
    sems.append(sem(x,n))

    plt.clf()
    plt.figure(figsize=(6,4))

    colors = sns.cubehelix_palette(len(fracs))
    pos = list(range(0,len(fracs)))
    pos[-1] = pos[-1]+0.5

    plt.bar(x = pos, height = fracs, yerr = sems,
            align='center',width=0.6,color=colors,alpha=.7,)

    plt.ylabel('Proportion of windows with\nconstraint Z>{0} (%)'.format(cutoff), fontsize=12)
    plt.xlabel('Percentage of coding bases (%)', fontsize=12)
    ticks = range(0,len(fracs))
    labels = list(range(0,110,10))+['Exonic\nonlly']
    ticks = range(0,len(labels))
    ticks = [i-0.5 for i in ticks]

    plt.xticks(ticks=ticks, labels=labels,)
    plt.yticks(ticks=np.arange(0, 0.16, .02),labels=[(int(i*100)) for i in np.arange(0, 0.16, .02)])

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 2c
def plt_prop_z4_cd_titration(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',
                       index_col='element_id')[['z','coding_prop']]
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')

    z_nc = df_z[df_z['coding_prop']==0]['z'].to_list()
    z_exonic = df_z_exonic['z'].to_list()

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
    plt.errorbar(y = fracs, x = range(0,len(fracs)), yerr=sems, 
                 color=color, ecolor=color, ls='', marker='o',elinewidth=2,alpha=0.7,)
    plt.ylabel('Proportion of windows with\nconstraint Z>{0} (%)'.format(cutoff), fontsize=12)
    plt.xlabel('Random draws of exonic windows (%)\ninto the non-coding genome', fontsize=12)

    plt.xticks(range(0,len(fracs)), np.arange(0,110,10))
    plt.yticks(ticks=np.arange(0.009, 0.0112, .0005),
               labels=['{:.2f}'.format((i*100)) for i in np.arange(0.009, 0.0112, .0005)])        

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 2d
def plt_z_ptcl_nc_cd(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('constraint_z_genome_1kb_exonic.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',
                       index_col='element_id')[['z','coding_prop']]
    df_z_exonic = pd.read_csv('fig_tables/constraint_z_genome_1kb_exonic.txt',sep='\t',index_col='element_id')

    z_nc = df_z[df_z['coding_prop']==0]['z'].to_list()
    z_exonic = df_z_exonic['z'].to_list()

    fracs = []
    for pctl in np.arange(10,110,10):
        z_pctl = np.percentile(z_exonic,pctl)
        pctl_nc = stats.percentileofscore(z_nc,z_pctl)
        fracs.append(pctl_nc)

    plt.clf()
    plt.figure(figsize=(6,4))
    
    color = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    color_ = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    plt.scatter(y = fracs, x = range(0,len(fracs)), color=color, alpha=0.7)

    plt.ylabel('Non-coding constraint Z percentile', fontsize=12)
    plt.xlabel('Exonic constraint Z percentile', fontsize=12)
    plt.xticks(range(0,len(fracs)), np.arange(10,110,10))
    plt.axvline(x = 8,color=color,linestyle='dashed',ymax=fracs[8]/100)
    plt.axhline(y = fracs[8],color=color_,linestyle='dashed',xmax=.9)
    plt.axvline(x = 4,color=color,linestyle='dashed',ymax=fracs[4]/100)
    plt.axhline(y = fracs[4],color=color_,linestyle='dashed',xmax=.5)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 3a
def plt_enrichment_gwas_exl_ccre(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')    
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    df_z = df_z[df_z['coding_prop']==0]
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    df_z['z_bin'] = pd.cut(df_z['z'], bins=cut_bins, labels=range((-1)*label,label))
    
    plt.clf()
    plt.figure(figsize=(5,4))

    ann1,ann2 = 'GWAS Catalog repl (ext)','ENCODE cCREs'
    for exl in [False, True]:
        if exl: 
            dfp=df_z[df_z[ann2]==0]
        else: dfp=df_z

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

        rightshift = 0.15 if exl else -0.15
        alpha = 0.6 if exl else 1.
        label = 'Exl. cCREs' if exl else 'Whole genome'
        fmt = '-' if exl else '--'
        color = sns.cubehelix_palette(start=.5, rot=-.5, )[-3]
        plt.errorbar(y = odds2plot, 
                     x = [i+0.5+rightshift for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o', ecolor=color, color=color,fmt=fmt,linestyle='', elinewidth=2, alpha=alpha,
                     label = label,
                    )

    plt.axhline(y = 1.0,color='#969696',linestyle='dashed')
    plt.rcParams['legend.title_fontsize'] = 12.
    plt.legend(loc='upper left', fontsize = 12., title='GWAS Catalog')
    plt.xticks(bins[1:],bins[1:])
    plt.xlabel('Constraint Z', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 3b
def plt_enrichment_gwas_ukb_exl_ccre(savefig):
    
    download_fig_table('UKBB_94traits_release1.traits')
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('ukb_fine-mapping_cs95.json')
    import json
    ukbb_traits = [line.strip().split('\t')[1] for line in 
                   open('fig_tables/UKBB_94traits_release1.traits').readlines()[1:]]
    trait2name = dict(line.strip().split('\t')[1:3] for line in 
                      open('fig_tables/UKBB_94traits_release1.traits').readlines()[1:])
    f = open('fig_tables/ukb_fine-mapping_cs95.json')
    d_ukbb = json.load(f)

    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t')
    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['pass_qc'])]
    all_elements = set(dfp['element_id'])
    z4_elements = set(dfp[dfp['z']>=4]['element_id'])
    odds2plot={}
    ci2plot={}
    pval2plot={}
    ci_val = {}
    for trait in ukbb_traits:
        hits = set(d_ukbb[trait])
        x1, n1 = len(all_elements.intersection(hits)),len(all_elements)
        x2, n2 = len(z4_elements.intersection(hits)),len(z4_elements)
        fc = (x2/n2)/(x1/n1) if x1>0 else 'inf'
        odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
        ci = ciOfOdds(x2,n2,x1,n1)
        trait = trait2name[trait]
        odds2plot[trait] = ci[0]
        ci2plot[trait] = [ci[0]-ci[1],ci[2]-ci[0]]
        pval2plot[trait] = pval
        ci_val[trait] = ci[1:]

    ann2 = 'ENCODE cCREs'
    dfp = dfp[~dfp[ann2]]
    all_elements_ = set(dfp['element_id'])
    z4_elements_ = set(dfp[dfp['z']>=4]['element_id'])
    odds2plot_={}
    ci2plot_={}
    pval2plot_={}
    ci_val_ = {}
    for trait in ukbb_traits:
        hits = set(d_ukbb[trait])
        x1, n1 = len(all_elements_.intersection(hits)),len(all_elements_)
        x2, n2 = len(z4_elements_.intersection(hits)),len(z4_elements_)
        fc = (x2/n2)/(x1/n1) if x1>0 else 'inf'
        odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
        ci = ciOfOdds(x2,n2,x1,n1)
        trait = trait2name[trait]
        odds2plot_[trait] = ci[0]
        ci2plot_[trait] = [ci[0]-ci[1],ci[2]-ci[0]]
        pval2plot_[trait] = pval
        ci_val_[trait] = ci[1:]
    l = sorted(ci_val.items(), key=lambda item: item[1][0], reverse=True) # sort by lower ci
    x = [i[0] for i in l if pval2plot[i[0]]<=0.05]
    
    y1 = [odds2plot[i] for i in x]
    y2 = [odds2plot_[i] for i in x]
    yerr1 = [ci2plot[i] for i in x]
    yerr2 = [ci2plot_[i] for i in x]
    yerr1 = [[ci2plot[i][0],0] for i in x] # only plot lower ci
    yerr2 = [[ci2plot_[i][0],0] for i in x]

    plt.clf()
    plt.figure(figsize=(4,len(y1)/2))
    offset=0.1
    color=sns.cubehelix_palette(start=.5, rot=-.5, )[-3]

    p1 = plt.errorbar(x = y1, y = [i+offset for i in range(0,len(x))], xerr=np.transpose(np.array(yerr1)), 
                      color=color, ecolor='#969696',linestyle='',marker='o',alpha=1.,elinewidth=2,)
    p2 = plt.errorbar(x = y2, y = [i-offset for i in range(0,len(x))], xerr=np.transpose(np.array(yerr2)), 
                      color=color, ecolor='#969696',linestyle='',marker='o',alpha=0.6,elinewidth=2,)
    plt.axvline(x = 1.0,color='#969696',linestyle='dashed')

    plt.xlabel('Enrichment of fine-mapped GWAS hits\nin constrained non-coding regions', fontsize=12)
    plt.ylabel('Disease/trait in UKB', fontsize=12)
    plt.yticks(range(0,len(x)),x)
    plt.legend([p1,p2],['Whole genome','Exl. cCREs'])
    plt.ylim(-1.5,len(x)+1.5)
    plt.xlim(-0.5,8)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 3c
def plt_gwas_ccre_z(savefig):
    
    download_fig_table('ccre_gwascatalog_z_1kb_nc_max.txt')
    df_ccre = pd.read_csv('fig_tables/ccre_gwascatalog_z_1kb_nc_max.txt', sep='\t')

    df_ccre['element_id'] = df_ccre['chrom']+'-'+df_ccre['start'].astype(str)+'-'+df_ccre['end'].astype(str)
    pctl = pd.qcut(df_ccre[['element_id','z']].set_index('element_id')['z'], 10,labels=range(1,11)[::-1]).to_dict()
    df_ccre['z_pctl'] = df_ccre['element_id'].map(pctl)

    df_ccre['element_len'] = df_ccre['end']-df_ccre['start']

    fracs = []
    sems = []
    p = sum(df_ccre['GWAS Catalog'])/sum(df_ccre['element_len'])
    d1 = df_ccre.groupby('z_pctl')['GWAS Catalog'].apply(sum).to_dict()
    d2 = df_ccre.groupby('z_pctl')['element_len'].apply(sum).to_dict()
    for k in d1:
        x,n = d1[k], d2[k],
        fracs.append(x/n)
        sems.append(sem(x,n))

    plt.clf()
    plt.figure(figsize=(6,4))
    x = np.arange(0, 10, 1)
    colors = sns.cubehelix_palette(15, start=.5, rot=-.5,)
    plt.bar(x, fracs[::-1], color = colors,width = 0.5, 
            align = 'center', edgecolor = None, ecolor='#252525')
    plt.ylabel('Incidence of GWAS hits (per bp)',fontsize=12.)
    plt.xlabel('cCRE constraint decile',fontsize=12.)
    plt.xticks(x, [str(i) for i in range(1,10+1)])

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        

# efig. 4a
def plt_z_vs_conservation_enh(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('genome_1kb.annot_dist2gene.txt')
    download_fig_table('genome_1kb.phastCons_mean.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t')
    d_dist2gene = dict([line.strip().split('\t')[3],float(line.strip().split('\t')[-1])] 
                       for line in open('fig_tables/genome_1kb.annot_dist2gene.txt').readlines())
    df_z['dist2gene'] = df_z['element_id'].map(d_dist2gene)
    d_phast = dict([line.strip().split('\t')[0],float(line.strip().split('\t')[-1])] for line in 
                   open('fig_tables/genome_1kb.phastCons_mean.txt').readlines()
                   if int(line.strip().split('\t')[2])>500)
    df_z['phastCons'] = df_z['element_id'].map(d_phast)

    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['dist2gene']>10*1000)]
    dfp = dfp.dropna(subset=['z','phastCons'])
    quant = 10
    dfp['z_quant'] = pd.qcut(dfp['z'],quant,labels=range(1,quant+1))
    dfp['phastCons_quant'] = pd.qcut(dfp['phastCons'],quant,labels=range(1,quant+1))

    ann = 'FANTOM enhancers'
    d_n = dfp.groupby(['z_quant','phastCons_quant'])[ann].apply(len).to_dict()
    d_x = dfp[dfp[ann]].groupby(['z_quant','phastCons_quant']).apply(len).to_dict()
    n0 = sum(d_n.values())
    x0 = sum(d_x.values())
    d_odds = {}
    d_pval = {}
    for k in d_n:
        x, n = d_x[k], d_n[k]
        odds, pval = stats.fisher_exact([[x, n-x], [x0, n0-x0]],alternative='two-sided')
        d_odds[k] = odds
        d_pval[k] = pval
    data = []
    labels = []
    for i in range(1,quant+1):
        data.append([d_odds[(j,i)] for j in range(1,quant+1)])
        l = []
        for j in range(1,quant+1):
            odds = '%.2f' % d_odds[(j,i)] if d_pval[(j,i)]<=0.05 else ''
            l.append(odds)
        labels.append(l)

    plt.clf()
    plt.figure(figsize=(7,5))
    with sns.axes_style('white'):
        ax = sns.heatmap(
            np.array(data[::-1]),
            square = True,
            vmin = 0., vmax = 2.,center = 1,
            cmap=sns.diverging_palette(220, 20, as_cmap=True),
            cbar=True,cbar_kws={'label': 'Enrichment of FANTOM enhancers'},
            linewidths = 0.2, linecolor= 'lightgrey',
            yticklabels = [str(i) for i in range(1,quant+1)][::-1],
            xticklabels = [str(i) for i in range(1,quant+1)],
            fmt = 's',
            annot = np.array(labels[::-1]),
                    )
        plt.yticks(rotation=0,fontsize=10.)
        plt.xlabel('Constraint Z decile',fontsize=12.)
        plt.ylabel('Conservation score decile',rotation=90,fontsize=12.)

        plt.tick_params(axis='both',top=False,right=False)
        
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 4b
def plt_z_vs_conservation_gwas(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('genome_1kb.annot_dist2gene.txt')
    download_fig_table('genome_1kb.phastCons_mean.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t')
    d_dist2gene = dict([line.strip().split('\t')[3],float(line.strip().split('\t')[-1])] 
                       for line in open('fig_tables/genome_1kb.annot_dist2gene.txt').readlines())
    df_z['dist2gene'] = df_z['element_id'].map(d_dist2gene)
    d_phast = dict([line.strip().split('\t')[0],float(line.strip().split('\t')[-1])] for line in 
                   open('fig_tables/genome_1kb.phastCons_mean.txt').readlines()
                   if int(line.strip().split('\t')[2])>500)
    df_z['phastCons'] = df_z['element_id'].map(d_phast)

    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['dist2gene']>10*1000)]
    dfp = dfp.dropna(subset=['z','phastCons'])
    quant = 10
    dfp['z_quant'] = pd.qcut(dfp['z'],quant,labels=range(1,quant+1))
    dfp['phastCons_quant'] = pd.qcut(dfp['phastCons'],quant,labels=range(1,quant+1))

    ann = 'GWAS Catalog'
    d_n = dfp.groupby(['z_quant','phastCons_quant'])[ann].apply(len).to_dict()
    d_x = dfp[dfp[ann]].groupby(['z_quant','phastCons_quant']).apply(len).to_dict()
    n0 = sum(d_n.values())
    x0 = sum(d_x.values())
    d_odds = {}
    d_pval = {}
    for k in d_n:
        x, n = d_x[k], d_n[k]
        odds, pval = stats.fisher_exact([[x, n-x], [x0, n0-x0]],alternative='two-sided')
        d_odds[k] = odds
        d_pval[k] = pval
    data = []
    labels = []
    for i in range(1,quant+1):
        data.append([d_odds[(j,i)] for j in range(1,quant+1)])
        l = []
        for j in range(1,quant+1):
            odds = '%.2f' % d_odds[(j,i)] if d_pval[(j,i)]<=0.05 else ''
            l.append(odds)
        labels.append(l)

    plt.clf()
    plt.figure(figsize=(7,5))
    with sns.axes_style('white'):
        ax = sns.heatmap(
            np.array(data[::-1]),
            square = True,
            vmin = 0., vmax = 2.,center = 1,
            cmap=sns.diverging_palette(220, 20, as_cmap=True),
            cbar=True,cbar_kws={'label': 'Enrichment of GWAS hits'},
            linewidths = 0.2, linecolor= 'lightgrey',
            yticklabels = [str(i) for i in range(1,quant+1)][::-1],
            xticklabels = [str(i) for i in range(1,quant+1)],
            fmt = 's',
            annot = np.array(labels[::-1]),
                    )
        plt.yticks(rotation=0,fontsize=10.)
        plt.xlabel('Constraint Z decile',fontsize=12.)
        plt.ylabel('Conservation score decile',rotation=90,fontsize=12.)

        plt.tick_params(axis='both',top=False,right=False)
        
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        

# efig. 5a
def plt_enrichment_re_10kb(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('genome_1kb.annot_dist2gene.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    d_dist2gene = dict([line.strip().split('\t')[3],float(line.strip().split('\t')[-1])] 
                       for line in open('fig_tables/genome_1kb.annot_dist2gene.txt').readlines())
    df_z['dist2gene'] = df_z['element_id'].map(d_dist2gene)

    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['dist2gene']>10*1000)]
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))

    plt.clf()
    plt.figure(figsize=(4,4))
    ann_color = { 
        'Super enhancers': sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
        'ENCODE cCREs': sns.cubehelix_palette(start=.5, rot=-.5, )[-3],
        'FANTOM enhancers':sns.cubehelix_palette(start=.5, rot=-.5, )[-4],
        'FANTOM lncRNAs': sns.cubehelix_palette(start=.5, rot=-.5, )[-5],
    }

    for ann in ['Super enhancers','ENCODE cCREs','FANTOM enhancers','FANTOM lncRNAs']:

        odds2plot = []
        ci2plot = []
        l1 = dfp[ann].dropna().to_list() # baseline
        x1,n1 = len([i for i in l1 if i]),len(l1)

        bins = sorted(list(set(dfp['z_bin'])))
        for z in bins:
            l2 = dfp[dfp['z_bin']==z][ann].dropna().to_list()
            x2,n2 = len([i for i in l2 if i]),len(l2)
            fc = (x2/n2)/(x1/n1) if x1>0 else 'inf'
            odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
            ci = ciOfOdds(x2,n2,x1,n1)
            odds2plot.append(ci[0])
            ci2plot.append([ci[0]-ci[1],ci[2]-ci[0]])
                
        color = ann_color[ann]
        plt.errorbar(y = odds2plot, x = [i+0.5 for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o',ecolor=color, color=color,fmt='--',elinewidth=2,alpha=0.7,label = ann,)

    plt.axhline(y = 1.0,color='#969696',linestyle='dashed')
    plt.legend(bbox_to_anchor=(1, 0.75),fontsize = 12.)
    plt.xticks(bins[1:],bins[1:],)
    plt.xlabel('Constraint Z', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 

# efig. 5b
def plt_enrichment_gwas_10kb(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    download_fig_table('genome_1kb.annot_dist2gene.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    d_dist2gene = dict([line.strip().split('\t')[3],float(line.strip().split('\t')[-1])] 
                       for line in open('fig_tables/genome_1kb.annot_dist2gene.txt').readlines())
    df_z['dist2gene'] = df_z['element_id'].map(d_dist2gene)

    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['dist2gene']>10*1000)]
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))

    plt.clf()
    plt.figure(figsize=(6, 4))
    
    ann_color = { 
        'GWAS Catalog': sns.cubehelix_palette(start=.5, rot=-.5, )[-5],
        'GWAS Catalog\n- repl (int)':  sns.cubehelix_palette(start=.5, rot=-.5, )[-4],
        'GWAS Catalog\n- repl (ext)':  sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
    }
    ann_label = {
        'GWAS Catalog': 'GWAS Catalog',
        'GWAS Catalog\n- repl (int)':  'Int repl',
        'GWAS Catalog\n- repl (ext)':  'Ext repl',
    }

    
    anns = ['GWAS Catalog','GWAS Catalog\n- repl (int)','GWAS Catalog\n- repl (ext)']
    for ann in anns:

        odds2plot = []
        ci2plot = []
        l1 = dfp[ann].dropna().to_list() # baseline
        x1,n1 = len([i for i in l1 if i]),len(l1)

        bins = sorted(list(set(dfp['z_bin'])))
        for z in bins:
            l2 = dfp[dfp['z_bin']==z][ann].dropna().to_list()
            x2,n2 = len([i for i in l2 if i]),len(l2)
            fc = (x2/n2)/(x1/n1) if x1>0 else 'inf'
            odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
            ci = ciOfOdds(x2,n2,x1,n1)
            odds2plot.append(ci[0])
            ci2plot.append([ci[0]-ci[1],ci[2]-ci[0]])

        rightshift = 0.15*anns.index(ann)
        color = ann_color[ann]
        label = ann_label[ann]
        plt.errorbar(y = odds2plot, x = [i+0.5+rightshift for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o',ecolor=color, color=color,fmt='--',linestyle='',elinewidth=2,alpha=0.7,
                     label = label,
                    )

    plt.axhline(y = 1.0,color='#969696',linestyle='dashed')
    plt.legend(loc='upper left', fontsize = 12.)
    plt.rcParams['legend.title_fontsize'] = 12.        
    plt.xticks(bins[1:],bins[1:],)
    plt.xlabel('Constraint Z', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 

# efig. 5c
def plt_enrichment_gwas_ukb_10kb(savefig):
    
    download_fig_table('UKBB_94traits_release1.traits')
    download_fig_table('ukb_fine-mapping_cs95.json')
    download_fig_table('genome_1kb.annot_dist2gene.txt')
    import json
    ukbb_traits = [line.strip().split('\t')[1] for line in 
                   open('fig_tables/UKBB_94traits_release1.traits').readlines()[1:]]
    trait2name = dict(line.strip().split('\t')[1:3] for line in 
                      open('fig_tables/UKBB_94traits_release1.traits').readlines()[1:])
    f = open('fig_tables/ukb_fine-mapping_cs95.json')
    d_ukbb = json.load(f)
    
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    d_dist2gene = dict([line.strip().split('\t')[3],float(line.strip().split('\t')[-1])] 
                       for line in open('fig_tables/genome_1kb.annot_dist2gene.txt').readlines())
    df_z['dist2gene'] = df_z['element_id'].map(d_dist2gene)

    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['dist2gene']>10*1000)]
    all_elements = set(dfp['element_id'])
    z4_elements = set(dfp[dfp['z']>=4]['element_id'])

    odds2plot={}
    ci2plot={}
    pval2plot={}
    ci_val = {}
    for trait in ukbb_traits:
        hits = set(d_ukbb[trait])
        x1, n1 = len(all_elements.intersection(hits)),len(all_elements)
        x2, n2 = len(z4_elements.intersection(hits)),len(z4_elements)
        fc = (x2/n2)/(x1/n1) if x1>0 else 'inf'
        odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
        ci = ciOfOdds(x2,n2,x1,n1)
        trait = trait2name[trait]
        odds2plot[trait] = ci[0]
        ci2plot[trait] = [ci[0]-ci[1],ci[2]-ci[0]]
        pval2plot[trait] = pval
        ci_val[trait] = ci[1:]

    l = sorted(ci_val.items(), key=lambda item: item[1][0], reverse=True) # sort by lower ci
    x = [i[0] for i in l if pval2plot[i[0]]<=0.05] 
    y = [odds2plot[i] for i in x]
    yerr = [[ci2plot[i][0],0] for i in x] # only plot lower ci for presentation
    
    plt.clf()
    plt.figure(figsize=(4, len(y)/3.5))

    plt.errorbar(y = x, x = y, xerr=np.transpose(np.array(yerr)), 
                 color=sns.cubehelix_palette(start=.5, rot=-.5, )[-3],ecolor='#969696',
                 linestyle='',marker='o',alpha=0.7,elinewidth=2,
                )

    plt.axvline(x = 1.0,color='#969696',linestyle='dashed')

    plt.xlabel('Enrichment of fine-mapped GWAS hits\nin constrained non-coding regions', fontsize=12)
    plt.ylabel('Disease/trait in UKB', fontsize=12)

    plt.ylim(-1.5,len(x)+1.5)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 6a
def plt_mave_enh_scores(savefig):     
    
    download_fig_table('mave_enh_scores.txt')
    df_mave = pd.read_csv("fig_tables/mave_enh_scores.txt",sep = "\t",index_col="element_name")
    scores = ["Constraint Z","Unfiltered Z","Orion","CDTS","gwRVIS","JARVIS","LINSIGHT"]

    plt.clf()
    fig, axes = plt.subplots(len(scores)+1, 1, figsize=(len(df_mave)*0.7,len(scores)*1.7),sharex="col")
    
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3]
    color3 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-4]
    colors = {"sig3":"#bdbdbd","sig2":color3,"sig1":color2}
    labels = {"sig3":"p>0.05","sig2":"1e-5<p<0.05","sig1":"p<1e-5"}

    x = range(0,len(df_mave))
    bottom_vals = [0]*len(x)
    for sig in ["sig1","sig2","sig3"]:
        vals = list(df_mave.sort_values(by="sig")[sig])
        axes[0].bar(x=x, height=vals, width = 0.7, color=colors[sig], edgecolor='white', alpha = .7,
                    bottom=list(bottom_vals), label = labels[sig])
        bottom_vals = [bottom_vals[i] + vals[i] for i in range(0,len(vals))]
    axes[0].set_ylabel("Percentage of\nmutations",fontsize = 12.)
    axes[0].legend(title="MAVE significance level",bbox_to_anchor=(1, 1), ncol=1,)
    plt.setp(axes[0], xticks=x, xticklabels=[i.split()[0] for i in list(df_mave.sort_values(by="sig").index)])
    axes[0].xaxis.set_tick_params(labelbottom=True, rotation=30)
    axes[0].set_title("Enhancers tested by MAVE",fontsize = 12.)
    
    for i in range(0,len(scores)):
        score = scores[i]
        ax = axes[i+1]  
        vals = list(df_mave.sort_values(by="sig")[score])
        ax.scatter(x=x, y=vals, color=color1, edgecolor='white', alpha = .7,s=50)
        ax.plot(x, vals, color=color1, linestyle='--', lw=1., alpha = .7)
        ax.set_ylabel(score,fontsize = 12.) 
        ax.xaxis.set_tick_params(labelbottom=False)
        dfp = df_mave[["sig",score]].dropna()
        print (len(dfp), score, 
               scipy.stats.pearsonr(dfp["sig"],dfp[score])[0], scipy.stats.spearmanr(dfp["sig"],dfp[score])[0])

    sns.despine(left=False, right=True, top=True, bottom=False)
    axes[0].tick_params(axis='both',top=False,right=False)
    fig.tight_layout() 

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 

# efig. 6b
def plt_mave_enh_scores_corr(savefig):     
    
    download_fig_table('mave_enh_scores.txt')
    df_mave = pd.read_csv("fig_tables/mave_enh_scores.txt",sep = "\t",index_col="element_name")
    scores = ["Constraint Z","Unfiltered Z","Orion","CDTS","gwRVIS","JARVIS","LINSIGHT"]
    cat = ['All',"Z scored","Exl. UC88"]

    D = {}
    for score in scores:
        dfp = df_mave[["sig",score]].dropna()
        pearson = scipy.stats.pearsonr(dfp["sig"],dfp[score])[0]
        spearman = scipy.stats.spearmanr(dfp["sig"],dfp[score])[0]
        D[score] = [spearman]
        dfp = df_mave[~df_mave["Constraint Z"].isna()][["sig",score]].dropna()
        pearson = scipy.stats.pearsonr(dfp["sig"],dfp[score])[0]
        spearman = scipy.stats.spearmanr(dfp["sig"],dfp[score])[0]
        D[score] += [spearman]    
        dfp = df_mave.drop(index=["UC88 enhancer"])[["sig",score]].dropna()
        pearson = scipy.stats.pearsonr(dfp["sig"],dfp[score])[0]
        spearman = scipy.stats.spearmanr(dfp["sig"],dfp[score])[0]
        D[score] += [spearman]

    plt.clf()
    fig, ax = plt.subplots(figsize=(len(cat)*2,len(scores)*0.6))

    with sns.axes_style("white"):
        ax = sns.heatmap(np.array([D[i] for i in scores]),
                         vmin = -1., vmax = 1.,center = 0,
                         cmap=sns.diverging_palette(220, 20, as_cmap=True),
                         cbar=True,cbar_kws={'label': 'Correlation'},
                         linewidths = 0.2, linecolor= "lightgrey",
                         yticklabels = scores,
                         xticklabels = False,
                         annot = True,
                        )

        plt.yticks(rotation=0,fontsize=12.)
        plt.xticks(rotation=0,fontsize=10.)
        plt.title("Enhancers tested by MAVE",fontsize=12.)
        plt.tight_layout()
        axT = ax.twiny()
        axT.set_xticks([1/6,0.5,1/6*5][:len(cat)])
        axT.set_xticklabels(cat,fontsize=12)
        sns.despine(left=True, right=True, top=True, bottom=True)
        ax.tick_params(axis='both',top=False,bottom=False,left=False,right=False)
        axT.tick_params(axis='both',top=False,bottom=False,left=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# efig. 7a
def plt_power_depl_1kb(savefig):     
    
    download_fig_table('powered_constraint_z4_by_depl.txt')
    df_depl = pd.read_csv('fig_tables/powered_constraint_z4_by_depl.txt', sep='\t')
    # window size
    ws = 1000
    dfp = df_depl[df_depl['Window_size'] == ws]

    step = 0.1
    depl_min, depl_max = 0.1, 0.5
    l = np.arange(depl_min, depl_max+step,step)
    if ws == 1000: l_e = [0.16, 0.36, 0.42]
    if ws == 100: l_e = [0.27, 0.57, 0.65]
    l = sorted(list(l)+l_e)
    l = [round(i,3) for i in l]

    label_sfx = {0.16: '(50th pctl exon)', 0.36: '(90th pctl exon)', 0.42: '(95th pctl exon)',
                0.27: '(50th pctl exon)', 0.57: '(90th pctl exon)', 0.65: '(95th pctl exon)'}
    colors = sns.cubehelix_palette(len(l), start=.5, rot=-.5,)

    plt.clf()
    plt.figure(figsize=(6,4))

    for depl in l:
        x = dfp[dfp['depl']==depl]['Sample_size'].to_list()
        y = dfp[dfp['depl']==depl]['Percent'].to_list()
        i = l.index(depl)
        color = colors[i]
        lw = 3
        label = '{0}%'.format(round(depl*100))
        ls = 'solid'
        if depl in label_sfx:
            label += '\n'+label_sfx[depl]
            ls = 'dotted'
        plt.plot(x, y, color=color, lw=lw, label=label, ls=ls )

    plt.legend(title='Depletion of variation', fontsize = 11.,bbox_to_anchor=(1, 1), )
    plt.yticks([0,.2,.4,.6,.8, 1.0], ['0','20','40','60','80','100'])
    plt.xscale('log')
    plt.xlabel('Sample size required',fontsize = 12.,)
    if ws == 1000: plt.ylabel('Percent of 1kb windows powered\nto detect constraint (Z>4)',fontsize = 12.,)
    if ws == 100: plt.ylabel('Percent of 100bp windows powered\nto detect constraint (Z>4)',fontsize = 12.,)

    n_v3 = 76156
    plt.axvline(x = n_v3, color='#969696', linestyle='dashed',lw=2.0)
    plt.axhline(y = .8, color='#fcbba1',)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

# efig. 7b
def plt_power_depl_100bp(savefig):     
    
    download_fig_table('powered_constraint_z4_by_depl.txt')
    df_depl = pd.read_csv('fig_tables/powered_constraint_z4_by_depl.txt', sep='\t')
    # window size
    ws = 100
    dfp = df_depl[df_depl['Window_size'] == ws]

    step = 0.1
    depl_min, depl_max = 0.1, 0.5
    l = np.arange(depl_min, depl_max+step,step)
    if ws == 1000: l_e = [0.16, 0.36, 0.42]
    if ws == 100: l_e = [0.27, 0.57, 0.65]
    l = sorted(list(l)+l_e)
    l = [round(i,3) for i in l]

    label_sfx = {0.16: '(50th pctl exon)', 0.36: '(90th pctl exon)', 0.42: '(95th pctl exon)',
                0.27: '(50th pctl exon)', 0.57: '(90th pctl exon)', 0.65: '(95th pctl exon)'}
    colors = sns.cubehelix_palette(len(l), start=.5, rot=-.5,)

    plt.clf()
    plt.figure(figsize=(6,4))

    for depl in l:
        x = dfp[dfp['depl']==depl]['Sample_size'].to_list()
        y = dfp[dfp['depl']==depl]['Percent'].to_list()
        i = l.index(depl)
        color = colors[i]
        lw = 3
        label = '{0}%'.format(round(depl*100))
        ls = 'solid'
        if depl in label_sfx:
            label += '\n'+label_sfx[depl]
            ls = 'dotted'
        plt.plot(x, y, color=color, lw=lw, label=label, ls=ls )

    plt.legend(title='Depletion of variation', fontsize = 11.,bbox_to_anchor=(1, 1), )
    plt.yticks([0,.2,.4,.6,.8, 1.0], ['0','20','40','60','80','100'])
    plt.xscale('log')
    plt.xlabel('Sample size required',fontsize = 12.,)
    if ws == 1000: plt.ylabel('Percent of 1kb windows powered\nto detect constraint (Z>4)',fontsize = 12.,)
    if ws == 100: plt.ylabel('Percent of 100bp windows powered\nto detect constraint (Z>4)',fontsize = 12.,)

    n_v3 = 76156
    plt.axvline(x = n_v3, color='#969696', linestyle='dashed',lw=2.0)
    plt.axhline(y = .8, color='#fcbba1',)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)


# efig. 8a
def plt_hist_freq_z_chrx(savefig):
    
    download_fig_table('constraint_z_genome_1kb_chrx.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb_chrx.annot.txt',sep='\t',index_col='element_id')
    z1 = df_z[df_z['coding_prop']==0]['z'].to_list()
    z2 = df_z[df_z['coding_prop']>0]['z'].to_list()
    
    plt.clf()
    fig, ax1 = plt.subplots(figsize=(6,4))
    
    color1, color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], sns.cubehelix_palette(8)[2]
    label1, label2 = 'Non-coding','Coding'
    hist_bins = np.arange(-10,11,1)

    ax1.hist(z1, bins=hist_bins, edgecolor=None, color = color1, alpha=0.7,density=False, label = label1)
    ax1.hist([0], edgecolor=None, color = color2, alpha=0.7,density=False, label = label2)
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax2.hist(z2, bins=hist_bins, edgecolor=None, color = color2,alpha=0.7,density=False, label= label2)

    ax1.legend(loc='upper left', fontsize = 12.)
    ax1.set_xlabel('Constraint Z',fontsize = 12.)   
    ax1.set_ylabel('Frequency',fontsize = 12.,  color=color1)
    ax1.tick_params(axis='y', labelcolor=color1)
    ax2.set_ylabel('Frequency',fontsize = 12.,  color=color2)
    ax2.tick_params(axis='y', labelcolor=color2)
    
    plt.axvline(x = 0.0,color='#525252',linestyle='dashed')
    plt.xlim(-10,10)
    plt.xticks(range(-10,11),range(-10,11))
    sns.despine(left=False, right=False, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=True)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        
# efig. 8b
def plt_enrichment_re_chrx(savefig):
    
    download_fig_table('constraint_z_genome_1kb_chrx.annot.txt')    
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb_chrx.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    dfp = df_z[ (df_z['coding_prop']==0)]
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))

    plt.clf()
    plt.figure(figsize=(4,4))
    ann_color = { 
        'Super enhancers': sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
        'ENCODE cCREs': sns.cubehelix_palette(start=.5, rot=-.5, )[-3],
        'FANTOM enhancers':sns.cubehelix_palette(start=.5, rot=-.5, )[-4],
        'FANTOM lncRNAs': sns.cubehelix_palette(start=.5, rot=-.5, )[-5],
    }

    for ann in ['Super enhancers','ENCODE cCREs','FANTOM enhancers','FANTOM lncRNAs']:

        odds2plot = []
        ci2plot = []
        l1 = dfp[ann].dropna().to_list() # baseline
        x1,n1 = len([i for i in l1 if i]),len(l1)

        bins = sorted(list(set(dfp['z_bin'])))
        for z in bins:
            l2 = dfp[dfp['z_bin']==z][ann].dropna().to_list()
            x2,n2 = len([i for i in l2 if i]),len(l2)
            ci = ciOfOdds(x2,n2,x1,n1)
            odds2plot.append(ci[0])
            ci2plot.append([ci[0]-ci[1],ci[2]-ci[0]])
                
        color = ann_color[ann]
        plt.errorbar(y = odds2plot, 
                     x = [i+0.5 for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o',
                     ecolor=color, 
                     color=color,
                     fmt='--',
                     elinewidth=2,
                     alpha=0.7,
                     label = ann,
                    )

    plt.axhline(y = 1.0,color='#969696',linestyle='dashed')
    plt.legend(bbox_to_anchor=(1, 0.75),fontsize = 12.)
    plt.xticks(bins[1:],bins[1:],)
    plt.xlabel('Constraint Z', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 



