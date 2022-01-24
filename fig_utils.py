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


# Fig. 1a
def plt_hist_freq_z(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    z1 = df_z[df_z['coding_prop']==0]['z'].to_list()
    z2 = df_z[df_z['coding_prop']>0]['z'].to_list()
    
    plt.clf()
    fig, ax1 = plt.subplots()
    
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


# Fig. 1b
def plt_aps_vs_z(savefig):
    
    download_fig_table('gnomad_v2.1_sv.sites.lft38.z_1kb_bootstrap.txt')
    fh = open('fig_tables/gnomad_v2.1_sv.sites.lft38.z_1kb_bootstrap.txt').readlines()[1:]
    m = [float(line.strip().split('\t')[1]) for line in fh]
    bins = [int(line.strip().split('\t')[0]) for line in fh]
    lci = [float(line.strip().split('\t')[1])-float(line.strip().split('\t')[2]) for line in fh]
    hci = [float(line.strip().split('\t')[3])-float(line.strip().split('\t')[1]) for line in fh]

    plt.clf()

    color = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    off = 0.5
    bins = range(-10,10)

    plt.errorbar(y = m, x = [i+off for i in bins], yerr=[lci,hci], 
                 ecolor=color, color=color, ls='', marker='o', elinewidth=2, alpha=0.7)
    plt.axhline(y = 0.0,color='#525252',linestyle='dashed')

    plt.ylabel('APS',fontsize = 12.)
    plt.xlabel('Constraint Z',fontsize = 12.)
    plt.xticks(range(-10,11), range(-10,11))
    plt.xlim(-10,10)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 2a
def plt_enrichment_re(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['pass_qc'])]
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
            out = [z,x1,n1,x2,n2,fc]
            odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
            out.append((odds, pval))
            ci = ciOfOdds(x2,n2,x1,n1)
            out.append(ci)
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


# Fig. 2b
def plt_enrichment_9p21(savefig):


    download_fig_table('constraint_z_genome_1kb.annot.txt')
    ann = '9p21'
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    dfp = df_z[df_z['pass_qc']]
    z_cd = dfp[dfp['coding_prop']>0]['z'].to_list()
    z_oth = dfp[dfp[ann]==0]['z'].to_list()

    dfp = dfp[dfp['coding_prop']==0]
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))
    z_ann = dfp[dfp[ann]>0]['z'].to_list()

    odds2plot = []
    ci2plot = []
    l1 = dfp[ann].dropna().to_list() # baseline
    x1,n1 = len([i for i in l1 if i]),len(l1)
    bins = sorted(list(set(dfp['z_bin'])))
    for z in bins:
        l2 = dfp[dfp['z_bin']==z][ann].dropna().to_list()
        x2,n2 = len([i for i in l2 if i]),len(l2)
        fc = (x2/n2)/(x1/n1) if x1>0 else 'inf'
        out = [z,x1,n1,x2,n2,fc]
        odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
        out.append((odds, pval))
        ci = ciOfOdds(x2,n2,x1,n1)
        out.append(ci)
        odds2plot.append(ci[0])
        ci2plot.append([ci[0]-ci[1],ci[2]-ci[0]])

    plt.clf()
    fig, ax1 = plt.subplots()
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color2 = sns.cubehelix_palette(8)[2]
    color3 = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]

    ax1 = sns.kdeplot(z_cd,shade=True,color = color2,alpha=0.7,label = 'Coding')
    ax1 = sns.kdeplot(z_ann,shade=True,color = color1,alpha=0.7,label = ann)
    ax1.legend(loc='upper left', fontsize = 12.)
    ax1.set_ylabel('Density',fontsize = 12.)
    ax1.set_xlabel('Constraint Z',fontsize = 12.)
    plt.xlim(-10,10)
    plt.xticks(range(-10,11),range(-10,11))

    ax2 = ax1.twinx() # instantiate a second axes that shares the same x-axis
    ax2.set_ylabel('Enrichment', color=color3,fontsize = 12.)
    ax2.tick_params(axis='y', labelcolor=color3)
    ax2.errorbar(y = odds2plot, x = [i+0.5 for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                 ecolor=color3, 
                 color=color3,
                 fmt='--',
                 marker='o',
                 alpha=0.7,
                )
    ax2.axhline(y = 1.0,color='#969696',linestyle='dashed')

    sns.despine(left=False, right=False, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=True)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 2c
def plt_9p21_locus(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    ann = '9p21'
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})

    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['pass_qc'])]
    dfp = dfp[dfp[ann]>0][['element_id','z','9p21','9p21 enhancers','9p21 gwasCatalog','ENCODE cCREs',]]

    dfp[['chrom','start','end']] = dfp['element_id'].str.rsplit('-',expand=True)
    dfp['start'] = dfp['start'].astype(int)
    dfp['end'] = dfp['end'].astype(int)

    for ann2 in ['9p21 gwasCatalog','9p21 enhancers']:
        l1 = dfp[dfp[ann2]]['z']
        l2 = dfp[~dfp[ann2]]['z']

    # fill in gaps
    l = dfp['start'].to_list()
    l2 = np.arange(l[0],l[-1]+1000,1000)
    l2 = [i for i in l2 if i not in l]
    d = {}
    for k in dfp.columns: d[k] = [0]*len(l2)
    d['start'] = l2
    dfp = pd.concat([dfp, pd.DataFrame.from_dict(d)], ignore_index=True).sort_values(by='start')
    ann2 = '9p21 enhancers'
    ann3 = '9p21 gwasCatalog'
    dfp['row_num'] = np.arange(len(dfp))
    if ann2:
        dfp['marker'] = np.where(dfp[ann2], dfp['z'], 0)
    if ann3:
        dfp['marker2'] = np.where(dfp[ann3], dfp['z']+0.5, 0)
        dfp['marker2'] = np.where( (dfp[ann3]) & (dfp['z']<0), dfp['marker2']-1, dfp['marker2'])

    plt.clf()
    fig = plt.figure(figsize=(len(dfp)/10,4))
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3]
    color3 = sns.cubehelix_palette(8)[2]
    alpha=1.

    plt.bar(x = dfp['row_num'], height = dfp['z'],align='center',width=0.8,color=color1,alpha=alpha,)
    if ann2:
        plt.bar(x = dfp['row_num'], height = dfp['marker'],align='center',width=0.8, color=color2,alpha=alpha,
                label='ENCODE cCREs')
    if ann3:
        plt.scatter(x = dfp[dfp[ann3]>0]['row_num'], y = dfp[dfp[ann3]>0]['marker2'], marker='x', 
                    color=color3, alpha=alpha,
                    label='CAD/T2D GWAS')
    plt.ylim(-10,10)
    plt.legend(loc='upper left', fontsize = 12.)

    plt.ylabel('Constraint Z', fontsize=12)
    plt.xlabel('9p21 Coordinates (Mb)', fontsize=12)
    plt.xticks(ticks=dfp['row_num'].to_list()[::10], 
               labels=[(i/1000000) for i in dfp['start'].to_list()][::10],rotation=45)


    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        

# Fig. 3a
def plt_enrichment_gwas(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')    
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t').rename(
        columns={'GWAS Catalog repl (int)':'GWAS Catalog\n- repl (int)',
                 'GWAS Catalog repl (ext)':'GWAS Catalog\n- repl (ext)'})
    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['pass_qc'])]
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
            out = [z,x1,n1,x2,n2,fc]
            odds, pval = stats.fisher_exact([[x2, n2-x2], [x1, n1-x1]],alternative='greater')
            out.append((odds, pval))
            ci = ciOfOdds(x2,n2,x1,n1)
            out.append(ci)
            odds2plot.append(ci[0])
            ci2plot.append([ci[0]-ci[1],ci[2]-ci[0]])

        rightshift = 0.15*anns.index(ann)
        color = ann_color[ann]
        label = ann_label[ann]
        plt.errorbar(y = odds2plot, 
                     x = [i+0.5+rightshift for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o',
                     ecolor=color, 
                     color=color,
                     fmt='--',
                     linestyle='',
                     elinewidth=2,
                     alpha=0.7,
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


# Fig. 3b
def plt_enrichment_gwas_ukb(savefig):

    download_fig_table('UKBB_94traits_release1.traits')
    download_fig_table('ukb_fine-mapping_cs95.json')
    download_fig_table('constraint_z_genome_1kb.annot.txt')
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

    l = sorted(ci_val.items(), key=lambda item: item[1][0], reverse=True) # sort by lower ci
    x = [i[0] for i in l if pval2plot[i[0]]<=0.05] 
    y = [odds2plot[i] for i in x]
    yerr = [[ci2plot[i][0],0] for i in x] # only plot lower ci for presentation
    
    plt.clf()
    plt.figure(figsize=(4, len(y)/3.5))

    plt.errorbar(y = x, 
                 x = y, 
                 xerr=np.transpose(np.array(yerr)), 
                 color=sns.cubehelix_palette(start=.5, rot=-.5, )[-3],
                 ecolor='#969696',
                 linestyle='',
                 marker='o',
                 alpha=0.7,
                 elinewidth=2,
                )

    plt.axvline(x = 1.0,color='#969696',linestyle='dashed')

    plt.xlabel('Enrichment of fine-mapped GWAS hits\nin constrained non-coding regions', fontsize=12)
    plt.ylabel('Disease/trait in UKB', fontsize=12)

    plt.ylim(-1.5,len(x)+1.5)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 4a
def plt_dd_cnv_z(savefig):

    download_fig_table('cnvDevDelay_z_1kb_nc_max.txt')
    df_cnv = pd.read_csv('fig_tables/cnvDevDelay_z_1kb_nc_max.txt', sep = '\t', 
                         header=None, names = ['chrom','start','end','z','group'])
    d = df_cnv.groupby('group')['z'].apply(list).to_dict()
    
    counts = []
    for group in ['DD_control','DD_case','DD_pathogenic','ClinVar_pathogenic']:
        l = d[group]
        counts.append([len([i for i in l if i>=4]), len(l)])
    counts[2][0]+=1 # [15,18] --> [16,18] to include 'chr3  195988732   197628732   3.991431383' (z=3.991431383)
    props = [i[0]/i[1] for i in counts]
    sems = [sem(i[0],i[1]) for i in counts]

    plt.clf()
    fig = plt.figure(figsize=(4,4))
    
    plt.bar(x = range(0,len(props)), height = props, yerr = sems,
            align='center',width=0.6,color=sns.cubehelix_palette(8, start=.5, rot=-.5,)[2],alpha=.7,)
    plt.ylabel('Proportion of constrained CNVs', fontsize=12)
    plt.xticks(ticks=range(0,len(props)), 
               labels=['DD\ncontrol','DD\ncase','DD\npathogenic', 'ClinVar'],
               fontsize=12,)
    plt.ylim(0,1)
    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 4b
def plt_dd_cnv_logit(savefig):
    
    def cnv_dd_logit(y, x):
        import statsmodels.api as sm
        X = sm.add_constant(x, prepend=False)
        result = sm.Logit(y, X).fit()
        return [dict(result.params),dict(result.bse)]

    download_fig_table('cnvDevDelay_z_1kb_nc_max.logit.txt')
    df_logit = pd.read_csv('fig_tables/cnvDevDelay_z_1kb_nc_max.logit.txt', sep = '\t')
    
    ft = ['Non-coding constraint', 'Gene constraint','Gene number','CNV size']
    d_or = {}
    d_se = {}
    dfp = df_logit.drop(columns=['CNV type']).drop_duplicates()
    d_or['All'],d_se['All'] = cnv_dd_logit(dfp['DD case'], dfp[ft])
    dfp = df_logit[df_logit['CNV type'] == 'loss'].drop_duplicates()
    d_or['Deletion'],d_se['Deletion'] = cnv_dd_logit(dfp['DD case'], dfp[ft])
    dfp = df_logit[df_logit['CNV type'] == 'gain'].drop_duplicates()
    d_or['Duplication'],d_se['Duplication'] = cnv_dd_logit(dfp['DD case'], dfp[ft])
    
    plt.clf()
    fig,ax = plt.subplots(1, figsize=(4,4))
    mt = ['All', 'Deletion','Duplication']
    d_color = {'All': sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3], 
               'Deletion': sns.cubehelix_palette(8, start=.5, rot=-.5,)[-5], 
               'Duplication': sns.cubehelix_palette(8, start=.5, rot=-.5,)[-7]}   
    offset = 0
    for m in mt[::-1]:
        ors = []
        cis = []
        colors = []
        yt = []
        for f in ft:
            ors.append(d_or[m][f])
            ci = d_se[m][f]*1.96
            cis.append(ci)
            colors.append(d_color[m])
        ors = ors[::-1]
        cis = cis[::-1]
        colors = colors[::-1]
        yt = yt[::-1]
        alpha=.7
        plt.scatter(ors, np.arange(0+offset, (len(ors)+offset)*4, 4)[:4], 
                    s = 60, color = colors, alpha=alpha, label=m)
        plt.errorbar(x = ors, y = np.arange(0+offset, (len(ors)+offset)*4, 4)[:4], xerr=cis,  
                     ecolor=colors, fmt='None',alpha=alpha)
        offset +=1

    plt.legend(loc='best', fontsize = 12.)
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[2], handles[1], handles[0]]
    labels = [labels[2], labels[1], labels[0]]
    ax.legend(handles,labels,fontsize = 12.,loc='best',title='CNV type')

    plt.xlabel('DD Case vs. control log(OR)',fontsize = 12.)
    plt.yticks([1,5,9,13], ft[::-1], fontsize = 12.)
    plt.axvline(x=0., color='#969696')
    plt.grid(linestyle='--',axis='x')

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 4c
def plt_cnv_ihh(savefig):

    download_fig_table('cnv_IHH_dup4_z_1kb.txt')
    dfp = pd.read_csv('fig_tables/cnv_IHH_dup4_z_1kb.txt',sep='\t')

    ann2 = 'IHH_enh'
    dfp['row_num'] = np.arange(len(dfp))
    dfp['marker'] = np.where(dfp[ann2] > 0., dfp['z'], 0)

    plt.clf()
    fig,ax = plt.subplots(1, figsize=(len(dfp)/10,4))

    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3]

    ax.bar(x = dfp['row_num'], height = dfp['z'],align='center',width=0.8,color=color1)
    ax.bar(x = dfp['row_num'], height = dfp['marker'],align='center',width=0.8, color=color2,
           label='Major IHH enhancer')

    position = -3.2
    spacing = 0.3
    for n in [2,3,1,4]:
        dup_range = dfp[dfp['IHH_dup{0}'.format(n)]>0]['row_num'].to_list()
        position -= spacing
        plt.hlines(y = position, xmin=min(dup_range), xmax=max(dup_range), color='#969696',linewidth=3)
    n=4
    dup_range = dfp[dfp['IHH_dup{0}'.format(n)]>0]['row_num'].to_list()
    plt.hlines(y = position, xmin=min(dup_range), xmax=max(dup_range),color='#969696',linewidth=3,
               label='Dup in syndactyly/craniosynostosis')

    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[1], handles[0]]
    labels = [labels[1], labels[0]]
    ax.legend(handles,labels,fontsize = 12.,loc='upper left')

    plt.ylabel('Constraint Z', fontsize=12)
    plt.xlabel('2q35 Coordinates (Mb)', fontsize=12)
    plt.xticks(ticks=dfp['row_num'].to_list()[::10], 
               labels=[(i/1000000) for i in dfp['start'].to_list()][::10],rotation=45)
    plt.ylim(-6,6)


    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)


    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5a
def plt_prop_roadmaplinks(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t')
    dfp = df_z[ (df_z['coding_prop']==0) & (df_z['pass_qc'])]
    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))

    d = dfp.groupby('z_bin')['RoadmapLinks'].apply(list).to_dict()
    fracs = []
    sems = []
    for i in d:
        x, n = len([j for j in d[i] if j]),len(d[i])
        fracs.append(x/n)
        sems.append(sem(x,n))

    plt.clf()
    plt.figure(figsize=(6, 4))

    colors = sns.cubehelix_palette(max(len(fracs),15), start=.5, rot=-.5,)
    x = np.arange(0, len(fracs), 1)
    plt.bar(x, fracs, 
            color = colors,
            yerr = sems,
            width = 0.5, 
            align = 'center', edgecolor = None, ecolor='#252525')
    plt.ylabel('Proportion of non-coding\nwindows linked to a gene',fontsize=12.)
    plt.xlabel('Constraint Z',fontsize=12.)
    plt.xticks([i+0.5 for i in x[:-1]], cut_bins[1:-1], fontsize = 10,)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5b
def plt_enh_geneset_z(savefig):

    download_fig_table('enh_gene_roadmaplinks.txt')
    df_ge = pd.read_csv('fig_tables/enh_gene_roadmaplinks.txt',sep='\t')
    anns = ['Haploinsufficient','MGI essential','OMIM dominant','LOEUF constrained','LOEUF unconstrained','Olfactory']
    LL = []
    for dist2gene in [0,100*1000]:
        dfp = df_ge[df_ge['enh_gene_distance']>=dist2gene]
        L = []
        for ann in anns:
            l = dfp[dfp[ann]]['enhancer_constraint_Z'].to_list()
            L.append(l)
        LL.append(L)
        g1 = (dfp['Haploinsufficient']) | (dfp['MGI essential']) | (dfp['OMIM dominant']) | (dfp['LOEUF constrained'])
        g2 = (dfp['LOEUF unconstrained']) | (dfp['Olfactory'])
        l1 = dfp[g1]['enhancer_constraint_Z'].to_list()
        l2 = dfp[(~g1) & g2]['enhancer_constraint_Z'].to_list()

    plt.clf()
    fig, axs = plt.subplots(1, 2, figsize=(7.5, 5), sharey=True)

    boxprops = dict(linestyle='-', linewidth=0, color='white')
    whiskerprops = dict(linestyle='--', linewidth=1.5, color='#525252')
    medianprops = dict(linestyle='-', linewidth=1., color='white')
    colors =[sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]]*4 + ['#bababa']*2

    box0 = axs[0].boxplot(LL[0][::-1], 
                          vert = False, notch=True, patch_artist=True,meanline=True,widths = 0.5,showfliers=False,
                          boxprops = boxprops,medianprops = medianprops,whiskerprops = whiskerprops,)
    box1 = axs[1].boxplot(LL[1][::-1], 
                          vert = False, notch=True, patch_artist=True,meanline=True,widths = 0.5,showfliers=False,
                          boxprops = boxprops,medianprops = medianprops,whiskerprops = whiskerprops,)
    for patch, color in zip(box0['boxes'], colors[::-1]): patch.set_facecolor(color)
    for patch, color in zip(box1['boxes'], colors[::-1]): patch.set_facecolor(color)

    axs[0].set_xlabel('Enhancer constraint Z',fontsize = 12.)
    axs[1].set_xlabel('Distal enhancer constraint Z',fontsize = 12.)
    plt.yticks(range(1,len(anns)+1), anns[::-1], fontsize=12.)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5c
def plt_geneset_loeuf0_enhz1(savefig):

    download_fig_table('geneset_erichment_loeuf0_enhz1.txt')
    df_david = pd.read_csv('fig_tables/geneset_erichment_loeuf0_enhz1.txt',sep='\t',header=None).dropna()
    df_david.columns = df_david.iloc[0]
    df_david = df_david.drop(df_david.index[0])
    df_david['FDR'] = df_david['FDR'].astype(float)
    df_david['-log10(FDR)'] = (-1)*np.log10(df_david['FDR'])
    df_david = df_david.sort_values(by='-log10(FDR)')

    d_pval = dict(zip(df_david['Cluster Name'],df_david['-log10(FDR)']))

    plt.clf()
    plt.figure(figsize=(3.5,5))

    plt.barh(y=range(0,len(d_pval)), width=d_pval.values(), height = 0.6, 
             color = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], edgecolor='white')

    plt.axvline(x=(-1)*np.log10(0.05), color='#fb6a4a', linestyle='--')  
    plt.xlabel('-log10(FDR)',fontsize = 12.)
    plt.yticks(range(0,len(d_pval)), d_pval.keys(),fontsize = 12.)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5d
def plt_phastcons_loeuf0_enhz1(savefig):
    
    download_fig_table('enh_gene_roadmaplinks.txt')
    df_ge = pd.read_csv('fig_tables/enh_gene_roadmaplinks.txt',sep='\t')
    dfp = df_ge[df_ge['LOEUF unconstrained']]
    l1 = dfp[dfp['enhancer_constraint_Z']>=2.23]['gene_phastCons'].to_list()
    l2 = dfp[dfp['enhancer_constraint_Z']<2.23]['gene_phastCons'].to_list()

    plt.clf()
    plt.figure(figsize=(2.5,5.))

    boxprops = dict(linestyle='-', linewidth=0, color='white')
    whiskerprops = dict(linestyle='--', linewidth=1.5, color='#525252')
    medianprops = dict(linestyle='-', linewidth=1., color='white')
    flierprops = dict(marker='o', markerfacecolor='r', markersize=12,linestyle='none', markeredgecolor='g')
    colors = [sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]] + ['#bababa']

    box = plt.boxplot([l1,l2],
                      notch=True, patch_artist=True,meanline=True,widths = 0.6,showfliers=False,
                      boxprops = boxprops,medianprops = medianprops,whiskerprops = whiskerprops,)
    for patch, color in zip(box['boxes'], colors):patch.set_facecolor(color)
    for patch, color in zip(box['fliers'], colors):patch.set_markeredgecolor(color)  

    plt.xticks([1,2,],['yes','no'],)
    plt.ylabel('Gene conservation (phastCons)',fontsize = 12.)
    plt.xlabel('Constrained enhancer',fontsize = 12.)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5e    
def plt_enhz_tissue_expr_corr(savefig):    
    
    def enhz_gene_expr_lm(y, x):
        import statsmodels.api as sm
        X = sm.add_constant(x, prepend=False)
        result = sm.OLS(y, X).fit()
        return [dict(result.params),dict(result.bse)]

    download_fig_table('enh_gene_roadmaplinks.tissue_expr_gtex.txt')
    df_expr = pd.read_csv('fig_tables/enh_gene_roadmaplinks.tissue_expr_gtex.txt',sep='\t')
    df_expr[['tissue_roadmaplinks','tissue_gtex']]
    df_expr.index = df_expr.tissue_roadmaplinks + '---' + df_expr.tissue_gtex
    
    d_coeff = {}
    d_ci = {}

    for tissue in set(df_expr.index):
        dfp = df_expr[df_expr.index==tissue]
        result = enhz_gene_expr_lm(dfp['expression'], dfp[['enh_z','LOEUF']])
        d_coeff[tissue.split('---')[0]] = result[0]['enh_z']
        d_ci[tissue.split('---')[0]] = result[1]['enh_z']*1.96

    plt.clf()
    plt.figure(figsize=(3.5,5))

    l = sorted(d_coeff.items(), key=lambda item: item[1], reverse=True)
    y = [i[0] for i in l]
    x = [d_coeff[i] for i in y]
    xerr = [d_ci[i] for i in y]

    plt.errorbar(y = y, x = x, xerr=np.transpose(np.array(xerr)), 
                 color=sns.cubehelix_palette(8, start=.5, rot=-.5,)[2],ecolor='#969696', 
                 linestyle='', marker='o', alpha=0.7,elinewidth=2,)

    plt.axvline(x = 0,color='#969696',linestyle='dashed')
    plt.xlabel('Linear regression beta', fontsize=12)
    plt.yticks(fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 





