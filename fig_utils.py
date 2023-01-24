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


# Fig. 1a
def plt_hist_freq_gnocchi(savefig):
    
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    dfp = df_z[(df_z['pass_qc'])]
    z1 = dfp[(dfp['coding_prop']==0)]['z'].to_list()
    z2 = dfp[(dfp['coding_prop']>0)]['z'].to_list()
    
    plt.clf()
    fig, ax1 = plt.subplots(figsize=(6,4))
    
    color1, color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], sns.cubehelix_palette(8)[2]
    label1, label2 = 'Non-coding','Coding'
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


# Fig. 1b
def plt_aps_vs_gnocchi(savefig):
    
    download_fig_table('gnomad_v2.1_sv.sites.lft38.z_1kb_bootstrap.txt')
    fh = open('fig_tables/gnomad_v2.1_sv.sites.lft38.z_1kb_bootstrap.txt').readlines()[1:]
    m = [float(line.strip().split('\t')[1]) for line in fh]
    bins = [int(line.strip().split('\t')[0]) for line in fh]
    lci = [float(line.strip().split('\t')[1])-float(line.strip().split('\t')[2]) for line in fh]
    hci = [float(line.strip().split('\t')[3])-float(line.strip().split('\t')[1]) for line in fh]
    
    plt.clf()
    plt.figure(figsize=(6,4))

    off = 0.5
    bins = range(-5,5)

    plt.errorbar(y = m, x = [i+off for i in bins], yerr=[lci,hci], 
                 ls='', marker='o', elinewidth=2,
                 color = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3], alpha=.7)
    
    plt.axhline(y = 0.0,color='#525252',linestyle='dashed')
    plt.ylabel('APS',fontsize = 12.)
    plt.xlabel('Gnocchi',fontsize = 12.)
    plt.xticks(range(-4,5), range(-4,5))
    plt.xlim(-5,5)
    
    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 2a
def plt_enrichment_re(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    dfp = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]

    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))

    plt.clf()
    plt.figure(figsize=(6,4))
    cmap = sns.cubehelix_palette(start=.5, rot=-.5, )
    cmap2 = sns.cubehelix_palette()
    ann_color = {   
        'ENCODE cCRE-PLS': cmap2[-2],
        'ENCODE cCRE-pELS': cmap2[-3],
        'ENCODE cCRE-dELS':cmap2[-4],
        'ENCODE CTCF-only':'#969696',
        'Super enhancers': cmap[-2],
        'FANTOM enhancers':cmap[-4],
    }

    for ann in ['ENCODE cCRE-PLS','ENCODE cCRE-pELS','ENCODE cCRE-dELS','ENCODE CTCF-only'] + ['Super enhancers','FANTOM enhancers']:

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
        plt.errorbar(y = odds2plot, 
                     x = [i+0.5 for i in bins], yerr=np.transpose(np.array(ci2plot)), 
                     marker='o',
                     ecolor=color, 
                     color=color,
                     fmt='--',
                     linestyle='',
                     elinewidth=2,
                     alpha=0.7,
                     label = ann,
                    )

    plt.axhline(y = 1.0,color='#969696',linestyle='dashed')
    plt.legend(loc='upper left', fontsize = 12.)
    plt.xticks(bins[1:],bins[1:],)
    plt.xlabel('Gnocchi', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 2b
def plt_enrichment_gwas(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    dfp = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]

    cut_bins = np.array([-10]+list(np.arange(-4,5,1))+[10])
    label = int((len(cut_bins)-1)/2)
    dfp['z_bin'] = pd.cut(dfp['z'], bins=cut_bins, labels=range((-1)*label,label))

    plt.clf()
    plt.figure(figsize=(6, 4))
    
    ann_color = { 
        'GWAS Catalog': sns.cubehelix_palette(start=.5, rot=-.5, )[-5],
        'GWAS Catalog repl (ext)':  sns.cubehelix_palette(start=.5, rot=-.5, )[-4],
        'GWAS fine-mapping':  sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
    }
    ann_label = {
        'GWAS Catalog': 'GWAS Catalog',
        'GWAS Catalog repl (ext)':  '    w/ an independent replication',
        'GWAS fine-mapping': 'GWAS fine-mapping'
    }

    anns = ['GWAS Catalog','GWAS Catalog repl (ext)','GWAS fine-mapping']
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

        rightshift = 0.15*(anns.index(ann)-1)
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
    plt.xlabel('Gnocchi', fontsize=12)
    plt.ylabel('Enrichment', fontsize=12)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 2c
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
    dfp = df_z[(df_z['pass_qc']) & (df_z['coding_prop']==0)]
    all_elements = set(dfp['element_id'])
    z4_elements = set(dfp[dfp['z']>=4.]['element_id'])

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
                 color=sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
                 ecolor='#969696',
                 linestyle='',
                 marker='o',
                 alpha=0.7,
                 elinewidth=2,
                )

    plt.axvline(x = 1.0,color='#969696',linestyle='dashed')
    plt.xlabel('Enrichment of fine-mapped variants\nin constrained non-coding regions', fontsize=12)
    plt.ylim(-1.5,len(x)+1.5)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 
        

# Fig. 2d
def plt_plg_gwas(savefig):
    
    download_fig_table('plg_cad_cs.txt')
    download_fig_table('wgEncodeGencodeBasicV32.exons.bed')
    download_fig_table('wgEncodeGencodeBasicV32.introns.bed')

    segs_e = ['-'.join(line.strip().split('\t')[:3]) for line in 
              open('fig_tables/wgEncodeGencodeBasicV32.exons.bed').readlines()
             if 'ENST00000308192' in line.strip().split('\t')[3]]
    segs_e = [i.split('-') for i in segs_e]
    segs_e = [(i[0],int(i[1]),int(i[2])) for i in segs_e]
    segs_e = sorted(segs_e, key = lambda x: x[1])
    segs_i = ['-'.join(line.strip().split('\t')[:3]) for line in 
              open('fig_tables/wgEncodeGencodeBasicV32.introns.bed').readlines()
             if 'ENST00000308192' in line.strip().split('\t')[3]]
    segs_i = [i.split('-') for i in segs_i]
    segs_i = [(i[0],int(i[1]),int(i[2])) for i in segs_i]
    segs_i = sorted(segs_i, key = lambda x: x[1])
    
    dfp = pd.read_csv('fig_tables/plg_cad_cs.txt',sep='\t')

    plt.clf()
    fig = plt.figure(figsize=(len(dfp)/10,3))
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3]
    color3 = sns.cubehelix_palette(8)[2]
    alpha=.7

    plt.bar(x = dfp[dfp['marker']==0]['row_num'], 
            height = dfp[dfp['marker']==0]['z'],align='center',width=0.8,color='#bdbdbd',alpha=0.7,)
    plt.bar(x = dfp[dfp['marker']>0]['row_num'], height = dfp[dfp['marker']>0]['marker'],align='center',width=0.8, color=color1,alpha=alpha,
            label='CAD GWAS fine-mapping')
    x = dfp[dfp['CAD_CS_n']>0]['row_num'].to_list()
    y = (dfp[dfp['CAD_CS_n']>0]['z']+0.2).to_list()
    s = dfp[dfp['CAD_CS_n']>0]['CAD_CS_n'].astype(int).astype(str).to_list()
    for i in range(0,len(x)):
        plt.text(x = x[i], y = y[i], s = s[i],color=color1, alpha=alpha,ha='center', va='bottom',fontsize=8)   

    position = -3.2
    position = min(dfp['z'])-1.5
    spacing = 1.

    seg_range = dfp[dfp['plg_introns']+dfp['plg_exons']>0]['row_num'].to_list()
    x = np.arange(min(seg_range)+1,max(seg_range),1.5)
    plt.scatter(x = x, y = [position]*len(x), marker='4', s=50, linewidths=0.5,color='#969696')
    plt.hlines(y = position, xmin=min(seg_range), xmax=max(seg_range), color='#969696',linewidth=2)
    
    for seg in segs_e:
        if seg[1]-min(dfp['start'])<0:continue
        plt.hlines(y = position, 
                   xmin=(seg[1]-min(dfp['start']))/1000, 
                   xmax=(seg[1]-min(dfp['start']))/1000+0.2, 
                   color='#737373',linewidth=8) 

    seg_range = dfp[dfp['plg_antisense_introns']>0]['row_num'].to_list()
    position -= spacing
    plt.hlines(y = position, xmin=min(seg_range), xmax=max(seg_range), color='#969696',linewidth=2,) 
    x = np.arange(min(seg_range)+1,max(seg_range),1.5)
    plt.scatter(x = x, y = [position]*len(x), marker='3', s=50, linewidths=0.5,color='#969696')

    plt.ylim(-6,12)
    plt.legend(loc='upper left', fontsize = 12.)
    plt.ylabel('Gnocchi', fontsize=12)
    plt.xlabel('Coordinates (Mb)', fontsize=12)
    plt.xticks(ticks=dfp['row_num'].to_list()[::10], 
               labels=[(i/1000000) for i in dfp['start'].to_list()][::10],rotation=45)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 3a-b
def plt_comparison_roc(pos,neg,dist2exon,savefig):
    
    # download_fig_table('comparisons_*.txt')
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')
    
    scores = ['Orion','CDTS','gwRVIS','DR','phastCons','phyloP','GERP']
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


# Fig. 3c-d
def plt_dominance_scores(pos,savefig):
    
    from dominance_analysis import Dominance
    # download_fig_table('comparisons_*.txt')
    download_fig_table('comparisons.tar.gz')
    os.chdir('fig_tables')
    os.system('tar -xvf comparisons.tar.gz')
    os.chdir('../')

    scores = ['z','Orion','CDTS','gwRVIS','DR','phastCons','phyloP','GERP']
    score_color = {
        'z':sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
        'Orion':'#33a02c','CDTS':'#fb9a99','DR':'#993404','gwRVIS':'#6baed6',   
        'phastCons':'#969696','phyloP':'#737373','GERP':'#bdbdbd'
    }
    score_label = dict([i,i] for i in scores)
    score_label['z'] = 'Gnocchi'

    df_pos = {'gwas_catalog': 'comparisons_gwas_catalog_repl', 
              'gwas_fine-mapping':'comparisons_gwas_fine-mapping_pip09', 
              'gwas_fine-mapping_hc':'comparisons_gwas_fine-mapping_pip09_hc', 
              'clinvar_plp_hgmd':'comparisons_likely_pathogenic_clinvar_hgmd',     
          }
    title_label = {'gwas_catalog':'GWAS Catalog','gwas_fine-mapping':'GWAS fine-mapping',
                   'gwas_fine-mapping_hc':'GWAS fine-mapping\n(high confidence)',
                   'clinvar_plp_hgmd':'Likely pathogenic'}
    df_neg = {
        'topmed_maf5':'comparisons_topmed_maf5.sampled.cov', 
        'topmed_mac1':'comparisons_topmed_mac1.sampled.cov', 
          }

    sampling = 10

    if pos == 'gwas':
        df_1 = pd.concat([pd.read_csv('fig_tables/comparisons/comparisons_gwas_catalog_repl.txt',sep='\t'),
                          pd.read_csv('fig_tables/comparisons/comparisons_gwas_fine-mapping_pip09.txt',sep='\t')
                         ]).drop_duplicates(subset=['locus'])
        neg = 'topmed_maf5'
    elif pos=='clinvar_plp_hgmd':
        df_1 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_pos[pos]),sep='\t')
        neg = 'topmed_mac1'

    df_0 = pd.read_csv('fig_tables/comparisons/{0}.txt'.format(df_neg[neg]),sep='\t')    
    df_1['group'] = 1
    df_0['group'] = 0

    if sampling: df_0 = df_0.sample(n=sampling*len(df_1), random_state=714)

    df_01 = pd.concat([df_1,df_0]).drop_duplicates(subset=['locus'])
    for score in scores:
        df_01['{0}_'.format(score)] = df_01[score]-np.nanmin(df_01[score])

    df_01 = df_01.dropna(subset=scores)
    df_y = df_01['group']
    df_x = df_01[[i+'_' for i in scores]]

    dft = pd.concat([df_x,df_y],axis=1)
    dominance_classification=Dominance(data=dft,target='group',objective=0,
                                       pseudo_r2='mcfadden',
                                       top_k=df_x.shape[1],
                                      )
    incr_variable_rsquare=dominance_classification.incremental_rsquare()
    dfp = dominance_classification.dominance_stats()

    dom = 'Percentage Relative Importance'

    d_score_type = {
        'Gnocchi':'Human lineage-specific constraint',
        'Orion':'Human lineage-specific constraint',
        'CDTS':'Human lineage-specific constraint',
        'gwRVIS':'Human lineage-specific constraint',
        'DR':'Human lineage-specific constraint',
        'phastCons':'Interspecies conservation',
        'phyloP':'Interspecies conservation',
        'GERP':'Interspecies conservation',
                   }
    c = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    c2 = sns.cubehelix_palette(start=.5, rot=-.5, )[-3]
    d_color = {'Human lineage-specific constraint':c,'Interspecies conservation':c2}
    d_yl = {'gwas':'GWAS','clinvar_plp_hgmd':'Likely pathogenic'}

    dfp['score'] = [score_label['_'.join(i.split('_')[:-1])] for i in dfp.index]
    dfp['score_type'] = dfp['score'].map(d_score_type)
    dfp['row_number'] = np.arange(len(dfp))

    plt.clf()
    plt.figure(figsize=(6,4))
    for score_type in ['Human lineage-specific constraint','Interspecies conservation']:
        plt.bar(x = dfp[dfp['score_type']==score_type]['row_number'], 
                height = dfp[dfp['score_type']==score_type][dom],
                align='center',width=0.8,color=d_color[score_type],label=score_type,alpha=0.7)

    plt.legend(fontsize=11)
    plt.ylabel('Relative contribution in\n {0} variant classification (%)'.format(d_yl[pos]), fontsize=12)
    plt.xticks(range(0,len(dfp)),list(dfp['score']), fontsize=12,rotation=90)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 4a
def plt_cnv_dd_gnocchi(savefig):

    download_fig_table('cnvDevDelay_z_1kb_nc_max.txt')
    df_cnv = pd.read_csv('fig_tables/cnvDevDelay_z_1kb_nc_max.txt', sep = '\t', 
                         header=None, names = ['chrom','start','end','z','group'])
    df_cnv['element_id'] = df_cnv['chrom']+'-'+df_cnv['start'].astype(str)+'-'+df_cnv['end'].astype(str)

    props_l = []
    sems_l = []
    for dfp in [df_cnv]:
        d = dfp.groupby('group')['z'].apply(list).to_dict()
        counts = []
        for group in ['DD_control','DD_case','DD_pathogenic','ClinVar_pathogenic']:
            if group not in d: 
                counts.append([0,1])
                continue
            l = d[group]
            counts.append([len([i for i in l if i>=4]), len(l)])
        
        props = [i[0]/i[1] for i in counts]
        sems = [sem(i[0],i[1]) for i in counts]
        

        props_l.append(props)
        sems_l.append(sems)

    plt.clf()
    fig = plt.figure(figsize=(4,4))

    plt.bar(x = range(0,len(props_l[0])),  height = props_l[0], yerr = sems_l[0],
            align='center',width=0.6,
            color=sns.cubehelix_palette(start=.5, rot=-.5, )[-2],alpha=.7,)

    plt.ylabel('Proportion of constrained CNVs (%)', fontsize=12)
    plt.xticks(ticks=range(0,len(props)), 
               labels=['DD\ncontrol','DD\ncase','DD\npathogenic', 'ClinVar'],
               fontsize=12,)
    plt.yticks(ticks=np.arange(0, 1.1, .2),labels=[(int(i*100)) for i in np.arange(0, 1.1, .2)])

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
    ft_ = ['Non-coding constraint', 'CNV size']

    d_or = {}
    d_se = {}
    dfp = df_logit.drop(columns=['CNV type']).drop_duplicates()
    d_or['All'],d_se['All'] = cnv_dd_logit(dfp['DD case'], dfp[ft])

    dfp = dfp[dfp['Gene number'] == 0]
    d_or['Non-coding'],d_se['Non-coding'] = cnv_dd_logit(dfp['DD case'], dfp[ft_])

    plt.clf()
    fig,ax = plt.subplots(1, figsize=(4,4))
    mt = ['All', 'Deletion','Duplication']
    mt = ['All', 'Non-coding']
    d_color = {'All': sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3], 
               'Deletion': sns.cubehelix_palette(8, start=.5, rot=-.5,)[-5], 
               'Duplication': sns.cubehelix_palette(8, start=.5, rot=-.5,)[-7],
               'All': sns.cubehelix_palette(start=.5, rot=-.5, )[-4],
               'Non-coding': sns.cubehelix_palette(start=.5, rot=-.5, )[-5],
               'All': sns.cubehelix_palette(start=.5, rot=-.5, )[-2],
               'Non-coding':sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
              }  

    offset = 0.5
    for m in mt[::-1]:
        ors = []
        cis = []
        colors = []
        yt = []
        for f in ft:
            if f not in d_or[m]: 
                ors.append(np.nan)
                cis.append(np.nan)
            else: 
                ors.append(d_or[m][f])
                ci = d_se[m][f]*1.96
                cis.append(ci)
            colors.append(d_color[m])
        ors = ors[::-1]
        cis = cis[::-1]
        colors = colors[::-1]
        yt = yt[::-1]
        alpha=0.7
        plt.scatter(ors, np.arange(0+offset, (len(ors)+offset)*4, 4)[:4], 
                    s = 60, color = colors, alpha=alpha, label=m)
        plt.errorbar(x = ors, y = np.arange(0+offset, (len(ors)+offset)*4, 4)[:4], xerr=cis,  
                     ecolor=colors, fmt='None',alpha=alpha)
        offset +=1

    plt.legend(loc='best', fontsize = 12.)
    handles,labels = ax.get_legend_handles_labels()
    handles = handles[::-1]
    labels = labels[::-1]
    ax.legend(handles,labels,fontsize = 12.,loc='best')

    plt.xlabel('DD case vs. control log(OR)',fontsize = 12.)
    plt.yticks([1,5,9,13], ['Non-coding constraint\n(Gnocchi)', 'Gene constraint\n(LOEUF)',
                            'Gene number','CNV size'][::-1], fontsize = 12.)
    plt.axvline(x=0., color='#969696')
    plt.grid(linestyle='--',axis='x')

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 4c
def plt_cnv_ihh(savefig):

    download_fig_table('cnv_IHH_dup4_z_1kb.txt')
    download_fig_table('genome_1kb.IHH_gene_coding.txt')
    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t')
    d = dict(zip(df_z['element_id'],df_z['z']))
    dfp = pd.read_csv('fig_tables/cnv_IHH_dup4_z_1kb.txt',sep='\t')
    dfp['z'] = dfp['element_id'].map(d)
    
    ihh = set([line.strip().split('\t')[-1] for line in 
               open('fig_tables/genome_1kb.IHH_gene_coding.txt').readlines()])
    dfp['IHH_gene'] = dfp['element_id'].isin(ihh)
    ann2 = 'IHH_enh'
    dfp['row_num'] = np.arange(len(dfp))
    dfp['marker'] = np.where(dfp[ann2] > 0., dfp['z'], 0)

    ann3 = 'IHH_gene'
    if ann3:
        dfp['row_num'] = np.arange(len(dfp))
        dfp['marker2'] = np.where(dfp[ann3] > 0., dfp['z'], 0)
    
    plt.clf()
    fig,ax = plt.subplots(1, figsize=(len(dfp)/10,4))
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-2]
    color1 = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color3 = sns.cubehelix_palette(8)[2]

    ax.bar(x = dfp[(dfp['marker'] + dfp['marker2']==0)]['row_num'], 
           height = dfp[ (dfp['marker'] + dfp['marker2'])==0]['z'],align='center',width=0.8,color=color2,alpha=.7)
    
    ax.bar(x = dfp[dfp['marker']>0]['row_num'], 
           height = dfp[dfp['marker']>0]['marker'],align='center',width=0.8, color=color1,alpha=.7,
           label='Major IHH enhancer',
          )
    if ann3: 
        ax.bar(x = dfp[dfp['marker2']>0]['row_num'], 
               height = dfp[dfp['marker2']>0]['marker2'],
               align='center',width=0.8, color=sns.cubehelix_palette(8)[2],
               label='IHH gene',alpha=0.7,
          )
        
    position = -3.
    spacing = 0.3
    for n in [2,3,1,4]:
        dup_range = dfp[dfp['IHH_dup{0}'.format(n)]>0]['row_num'].to_list()
        position -= spacing
        plt.hlines(y = position, xmin=min(dup_range), xmax=max(dup_range), color='#969696',linewidth=2)
    n=4
    dup_range = dfp[dfp['IHH_dup{0}'.format(n)]>0]['row_num'].to_list()
    plt.hlines(y = position, xmin=min(dup_range), xmax=max(dup_range),color='#969696',linewidth=2,
               label='Dup in syndactyly/craniosynostosis')

    # plt.legend(loc='upper right',fontsize = 12.)
    plt.ylabel('Gnocchi', fontsize=12)
    plt.xlabel('2q35 Coordinates (Mb)', fontsize=12)
    plt.xticks(ticks=dfp['row_num'].to_list()[::10], 
               labels=[(i/1000000) for i in dfp['start'].to_list()][::10],rotation=45)
    plt.ylim(-5,10)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 4d
def plt_cnv_recurrent(savefig):

    download_fig_table('cnv_dd_recurrent.txt')
    dfp = pd.read_csv('fig_tables/cnv_dd_recurrent.txt',sep='\t')

    l = dfp['start'].to_list()
    l2 = np.arange(l[0],l[-1]+1000,1000)
    l2 = [i for i in l2 if i not in l]
    d = {}
    for k in dfp.columns: d[k] = [0]*len(l2)
    d['start'] = l2
    dfp = pd.concat([dfp, pd.DataFrame.from_dict(d)], ignore_index=True).sort_values(by='start')
    dfp['row_num'] = np.arange(len(dfp))

    dfp['seg_cases'] = dfp[[i for i in dfp.columns if 'seg_case_' in i]].sum(axis=1)
    dfp['seg_ctrls'] = dfp[[i for i in dfp.columns if 'seg_ctrl_' in i]].sum(axis=1)

    plt.clf()
    fig = plt.figure(figsize=(len(dfp)/25,4))
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color2 = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[-3]
    color3 = sns.cubehelix_palette(8)[2]
    alpha=.7

    plt.bar(x = dfp[dfp['seg_cases']<12]['row_num'], 
            height = dfp[dfp['seg_cases']<12]['z'],align='center',width=0.8,color=color1,alpha=alpha,
           linewidth=0)
    plt.bar(x = dfp[dfp['seg_cases']==12]['row_num'], 
            height = dfp[dfp['seg_cases']==12]['z'],align='center',width=0.8,color=color2,alpha=alpha,
           linewidth=0,label='Potential critical region')

    position = -3.2
    position = min(dfp['z'])-0.05
    spacing = 0.3
    for idx in range(0,12):
        seg_range = dfp[dfp['seg_case_{0}'.format(idx+1)]>0]['row_num'].to_list()
        position -= spacing
        plt.hlines(y = position, xmin=min(seg_range), xmax=max(seg_range), 
                   color=sns.cubehelix_palette(15)[0],
                   linewidth=2)  
    pos1 = position

    for idx in range(0,2):
        seg_range = dfp[dfp['seg_ctrl_{0}'.format(idx+1)]>0]['row_num'].to_list()
        position -= spacing
        plt.hlines(y = position, xmin=min(seg_range), xmax=max(seg_range), color='#bdbdbd',linewidth=2,)
    pos2 = position


    idx=12-1
    seg_range = dfp[dfp['seg_case_{0}'.format(idx+1)]>0]['row_num'].to_list()
    plt.hlines(y = pos1, xmin=min(seg_range), xmax=max(seg_range), 
               color=sns.cubehelix_palette(15)[0],
               linewidth=2,
               label='Deletion in DD case') 
    idx=2-1
    seg_range = dfp[dfp['seg_ctrl_{0}'.format(idx+1)]>0]['row_num'].to_list()
    plt.hlines(y = pos2, xmin=min(seg_range), xmax=max(seg_range), color='#bdbdbd',linewidth=2,
               label='Deletion in control') 

    plt.ylim(-10,10)
    # plt.legend(loc='upper left', fontsize = 12.)

    plt.ylabel('Gnocchi', fontsize=12)
    plt.xlabel('Coordinates (Mb)', fontsize=12)
    plt.xticks(ticks=dfp['row_num'].to_list()[::20], 
               labels=[(i/1000000) for i in dfp['start'].to_list()][::20],rotation=45)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 

def plt_cnv_recurrent_gnocchi(savefig):
    
    download_fig_table('cnv_dd_recurrent.txt')
    dfp = pd.read_csv('fig_tables/cnv_dd_recurrent.txt',sep='\t')
    dfp['seg_cases'] = dfp[[i for i in dfp.columns if 'seg_case_' in i]].sum(axis=1)
    dfp['seg_ctrls'] = dfp[[i for i in dfp.columns if 'seg_ctrl_' in i]].sum(axis=1)
    
    plt.clf()
    fig,ax = plt.subplots(1, figsize=(4.,4))
    color1 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    color2 = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]

    l1 = list(dfp[ (dfp['seg_cases'] ==12)]['z'])
    l2 = list(dfp[ (dfp['seg_cases'] > 0) & (dfp['seg_cases'] < 12)]['z'])
    print (len(l1),len(l2),np.median(l1),np.median(l2),u(l1,l2,alternative='greater').pvalue)
    sns.kdeplot(l2,shade=True,alpha=0.7,
                color=color1,
                label = 'Case, N<12')
    sns.kdeplot(l1,shade=True,alpha=0.7,
                color=color2,
                label = 'Case, N=12')
    l = list(dfp[ (dfp['seg_ctrls'] > 0)]['z'])
    sns.kdeplot(l,shade=False,alpha=1.,color='#bdbdbd',label = 'Control',lw=2)

    handles,labels = ax.get_legend_handles_labels()
    leg = ax.legend(handles[:-1][::-1]+handles[-1:],labels[:-1][::-1]+labels[-1:],
                    loc='upper right',
             )
    leg._legend_box.align = 'left'
    plt.xlabel('Gnocchi', fontsize=12)
    plt.ylabel('Density', fontsize=12)
    plt.xlim(-6,10)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)
    
    if savefig:
        plt.savefig(savefig, bbox_inches='tight')


# Fig. 5a
def plt_prop_roadmaplinks(savefig):

    download_fig_table('constraint_z_genome_1kb.annot.txt')
    df_z = pd.read_csv('fig_tables/constraint_z_genome_1kb.annot.txt',sep='\t',index_col='element_id')
    
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
    plt.figure(figsize=(4.5, 4.5))

    colors = sns.cubehelix_palette(max(len(fracs),15), start=.5, rot=-.5,)
    x = np.arange(0, len(fracs), 1)
    plt.bar(x, [i*100 for i in fracs], 
            color = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2],alpha=0.7,
            yerr = [i*100 for i in sems],
            width = 0.5, 
            align = 'center', edgecolor = None, ecolor='#252525')
    plt.ylabel('Proportion of non-coding regions\nas enhancers of specific genes (%)',fontsize=12.)
    plt.xlabel('Gnocchi',fontsize=12.)
    plt.xticks([i+0.5 for i in x[:-1]], cut_bins[1:-1], fontsize = 10,)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5b
def plt_enh_geneset_gnocchi(savefig):

    download_fig_table('enh_gene_roadmaplinks.txt')
    df_ge = pd.read_csv('fig_tables/enh_gene_roadmaplinks.txt',sep='\t')
    anns = ['Haploinsufficient','MGI essential','OMIM dominant','LOEUF constrained',
            'Olfactory','LOEUF unconstrained','LOEUF underpowered']
    L = []
    for ann in anns:
        l = df_ge[df_ge[ann]]['enhancer_constraint_Z'].to_list()
        L.append(l)

    plt.clf()
    plt.figure(figsize=(3.5,4.5))

    boxprops = dict(linestyle='-', linewidth=0, color='white',alpha=0.7)
    whiskerprops = dict(linestyle='--', linewidth=1.5, color='#525252',alpha=0.7)
    medianprops = dict(linestyle='-', linewidth=2., color='white',alpha=0.7)
    colors =[sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]]*4 + ['#bababa']*2 + ['white']
    ecolors =[sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]]*4 + ['#bababa']*3

    box = plt.boxplot(L[::-1], 
                      vert = False, notch=True, patch_artist=True,meanline=True,widths = 0.56,showfliers=False,
                      boxprops = boxprops,medianprops = medianprops,whiskerprops = whiskerprops)
    for patch, ecolor in zip(box['boxes'], ecolors[::-1]): patch.set(color=ecolor, linewidth=1.5)
    for patch, color in zip(box['boxes'], colors[::-1]): patch.set_facecolor(color)

    plt.xlabel('Enhancer Gnocchi',fontsize = 12.)
    plt.yticks(range(1,len(anns)+1), [i.strip(' 2') for i in anns[::-1]], fontsize=12.)

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5c
def plt_enh_gnocchi_loeuf_roc(savefig):

    download_fig_table('enhz_loeuf_pred.txt')
    df_ge = pd.read_csv('fig_tables/enhz_loeuf_pred.txt',sep='\t')
    df_test = df_ge[df_ge['train_test'] == 'test']
    df_train = df_ge[df_ge['train_test'] == 'train']

    ft = ['LOEUF']
    df_y = df_train['constrained']
    df_x = df_train[ft]
    logit = sm.Logit(df_y, sm.add_constant(df_x), has_constant='add').fit_regularized(disp=0)
    ft = ['LOEUF','enhancer_constraint_Z']
    df_x = df_train[ft]
    logit2 = sm.Logit(df_y, sm.add_constant(df_x), has_constant='add').fit_regularized(disp=0)

    df_test['pred1'] = logit.predict(sm.add_constant(df_test['LOEUF'], has_constant='add'))
    df_test['pred2'] = logit2.predict(sm.add_constant(df_test[['LOEUF','enhancer_constraint_Z']], has_constant='add'))

    dft = df_test[df_test['LOEUF_underpowered']]

    color1 = sns.cubehelix_palette(start=.5, rot=-.5, )[-2]
    color2 = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2]
    cmap = sns.cubehelix_palette(start=.5, rot=-.5, )
    score_color = {'pred1':color2,'pred2':color1}
    score_label = {'pred1':'LOEUF','pred2':'LOEUF+Enhancer\nGnocchi'}
    
    scores = ['pred1','pred2']
    d_auc = {}
    d_auc_ = {}
    idx = 0

    plt.clf()
    fig,ax = plt.subplots(1, figsize=(4.5,4.5))
    for score in scores:      
        y_true = dft['constrained']
        y_probas = dft[score]
        fpr, tpr, _ = roc_curve(y_true,  y_probas)
        auc = roc_auc_score(y_true, y_probas)
        plt.plot(fpr,tpr,color=score_color[score],alpha=.7, label='{0} ({1})'.format(score_label[score],round(auc,3)))
        d_auc[idx] = auc
        d_auc_[score] = auc
        idx +=1 
    plt.plot([0, 1], [0, 1], linestyle='--',color='#bdbdbd')

    sorted_auc = [i[0] for i in sorted(d_auc.items(), key=lambda item: item[1])][::-1]
    handles,labels = ax.get_legend_handles_labels()
    handles = [handles[i] for i in sorted_auc]
    labels = [labels[i] for i in sorted_auc]
    ax.legend(handles,labels,loc='lower right')

    plt.ylabel('True positive rate', fontsize=12)
    plt.xlabel('False positive rate', fontsize=12)
    plt.title("Classification of constrained genes\n(LOEUF underpowered)")

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 


# Fig. 5d
def plt_enh_gnocchi_tissue_expr_corr(savefig):    
    
    def enhz_gene_expr_lm(y, x):
        import statsmodels.api as sm
        x = x.apply(scipy.stats.zscore)
        X = sm.add_constant(x, prepend=False,has_constant='add')
        result = sm.OLS(y, X).fit()
        return [dict(result.params),dict(result.bse)]

    download_fig_table('enh_gene_roadmaplinks.tissue_expr_gtex.txt')
    df_expr = pd.read_csv('fig_tables/enh_gene_roadmaplinks.tissue_expr_gtex.txt',sep='\t')
    df_expr = df_expr[df_expr['expression']>=np.log2(1+1)]
    df_expr.index = df_expr.tissue_roadmaplinks + '---' + df_expr.tissue_gtex
    
    d_coeff = {}
    d_ci = {}

    for tissue in set(df_expr.index):
        dfp = df_expr[df_expr.index==tissue]
        result = enhz_gene_expr_lm(dfp[['expression']], dfp[['enh_z','LOEUF']])
        d_coeff[tissue.split('---')[0]] = result[0]['enh_z']
        d_ci[tissue.split('---')[0]] = result[1]['enh_z']*1.96

    plt.clf()
    plt.figure(figsize=(3.5,4.5))

    l = sorted(d_coeff.items(), key=lambda item: item[1], reverse=True)
    y = [i[0] for i in l]
    x = [d_coeff[i] for i in y]
    xerr = [d_ci[i] for i in y]
    plt.errorbar(y = y, x = x, xerr=np.transpose(np.array(xerr)), 
                 ecolor='#bdbdbd', 
                 fmt='None',
                 alpha=0.7,elinewidth=2,)
    plt.scatter(y = y, x = x, 
                color = sns.cubehelix_palette(8, start=.5, rot=-.5,)[2], alpha=0.7,
               s = 60, 
               )
        
    plt.axvline(x = 0,color='#969696',linestyle='dashed')
    plt.xlabel('Gene expression ~ \nEnhancer Gnocchi', fontsize=12) ### paper
    plt.ylabel('Tissue', fontsize=12)
    plt.title("Prediction of gene expression\n")

    sns.despine(left=False, right=True, top=True, bottom=False)
    plt.tick_params(axis='both',top=False,right=False)

    if savefig:
        plt.savefig(savefig, bbox_inches='tight') 




