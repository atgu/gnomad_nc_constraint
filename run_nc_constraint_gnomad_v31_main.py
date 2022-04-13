#!/usr/bin/env python
# coding: utf-8

import argparse
import pandas as pd
import csv
import os
import sys
import numpy as np
import scipy
from scipy import stats
import math
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score
from sklearn.model_selection import KFold
import pickle

input_bucket = 'gs://gnomad-nc-constraint-v31-paper'
cwd = os.getcwd()
sys.path.append(cwd)
os.system('gsutil cp {0}/misc/*.py {1}'.format(input_bucket,cwd))
from generic import *
from constraint_basics import *
from nc_constraint_utils import *
import hail as hl
hl.init()

def main(args):

    output_bucket = args.output_bucket
    output_dir = args.output_dir
    # skip_hail = args.skip_hail
    os.mkdir('{0}/tmp'.format(output_dir))

    ### Prefilter context ht and genome ht
    context_ht = hl.read_table('{0}/context_prepared.ht'.format(input_bucket))
    genome_ht = hl.read_table('{0}/genome_prepared.ht'.format(input_bucket))
    filter_to_autosomes_par(remove_coverage_outliers(context_ht)).write('{0}/context_prefiltered.ht'.format(output_bucket))
    filter_to_autosomes_par(remove_coverage_outliers(genome_ht)).write('{0}/genome_prefiltered.ht'.format(output_bucket))


    ### Compute mu for each context and methyl level from downsampled data

    # count observed and possible variants by context and methyl
    context_ht = hl.read_table('{0}/context_downsampled_1000.ht'.format(input_bucket))
    genome_ht = hl.read_table('{0}/genome_downsampled_1000.ht'.format(input_bucket))

    observed_ht = count_variants(genome_ht, count_downsamplings=['global', 'nfe', 'afr'],
                                 additional_grouping=('methyl_level',), omit_methylation=True)
    # observed_ht.write('{0}/observed_counts_by_context_methyl_downsampled_1000.ht'.format(output_bucket))

    grouping = hl.struct(context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, 
                         methylation_level=context_ht.methyl_level)
    output = {'variant_count': hl.agg.count()}
    possible_ht = context_ht.group_by(**grouping).aggregate(**output)
    # possible_ht.export('{0}/possible_counts_by_context_methyl_downsampled_1000.txt'.format(output_bucket))

    # compute mu
    observed_ht = observed_ht.annotate(possible_variants = possible_ht[observed_ht.key].variant_count)
    total_bases = observed_ht.aggregate(hl.agg.sum(observed_ht.possible_variants)) // 3
    total_mu = 1.2e-08

    correction_factors = observed_ht.aggregate(total_mu / (hl.agg.array_sum(observed_ht.downsampling_counts_global) / total_bases))
    observed_ht = annotate_variant_types(observed_ht.annotate(
        downsamplings_frac_observed=observed_ht.downsampling_counts_global / observed_ht.possible_variants,
        downsamplings_mu_snp=hl.literal(correction_factors) * observed_ht.downsampling_counts_global / observed_ht.possible_variants
        ))

    downsamplings = list(map(lambda x: x[1], get_downsamplings(observed_ht)))
    index_1kg = downsamplings.index(1000)

    observed_ht = observed_ht.annotate(
        observed_1kg=observed_ht.downsampling_counts_global[index_1kg],
        proportion_observed_1kg=observed_ht.downsampling_counts_global[index_1kg] / observed_ht.possible_variants,
        mu=observed_ht.downsamplings_mu_snp[index_1kg]
        )

    observed_ht.select(
        'transition','cpg','variant_type','variant_type_model',
        'possible_variants','observed_1kg','proportion_observed_1kg',
        'mu',
        ).export('{0}/mu_by_context_methyl_downsampled_1000.txt'.format(output_bucket))


    ### Compute proportion of possible variants observed for each context and methyl level from the whole dataset

    context_ht = hl.read_table('{0}/context_prefiltered.ht'.format(output_bucket))
    genome_ht = hl.read_table('{0}/genome_prefiltered.ht'.format(output_bucket))
    context_ht = filter_black_regions(context_ht.filter((context_ht.coverage_mean >= 30) & (context_ht.coverage_mean <= 32)))  
    genome_ht = genome_ht.semi_join(context_ht)

    # select & join
    context_ht = context_ht.select('context', 'ref', 'alt', 'methyl_level')
    genome_ht = genome_ht.select('context', 'ref', 'alt', 'methyl_level','freq', 'pass_filters')
    genome_join = genome_ht[context_ht.key]

    # filter for rare & pass qc variants
    af_cutoff = 0.001
    context_ht = context_ht.filter(hl.is_missing(genome_join) | ((genome_join.freq[0].AF <= af_cutoff) & genome_join.pass_filters))
    genome_ht = genome_ht.filter((genome_ht.freq[0].AF <= af_cutoff) & genome_ht.pass_filters)

    # count observed and possible variants by context and methyl
    grouping = hl.struct(context=genome_ht.context, ref=genome_ht.ref, alt=genome_ht.alt, methylation_level=genome_ht.methyl_level)
    output = {'variant_count': hl.agg.count()}
    observed_ht = genome_ht.group_by(**grouping).aggregate(**output)
    observed_ht.export('{0}/observed_counts_by_context_methyl.txt'.format(output_bucket))

    grouping = hl.struct(context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, methylation_level=context_ht.methyl_level)
    output = {'variant_count': hl.agg.count()}
    possible_ht = context_ht.group_by(**grouping).aggregate(**output)
    possible_ht.export('{0}/possible_counts_by_context_methyl.txt'.format(output_bucket))


    ### Fit proportion_observed ~ mu to obtain the context-specific mutabilities

    os.system('gsutil cp {0}/mu_by_context_methyl_downsampled_1000.txt {1}'.format(output_bucket,output_dir))
    os.system('gsutil cp {0}/observed_counts_by_context_methyl.txt {1}'.format(output_bucket,output_dir))
    os.system('gsutil cp {0}/possible_counts_by_context_methyl.txt {1}'.format(output_bucket,output_dir))
    def sem(x,n):
        import numpy as np
        p=float(x)/n
        return np.sqrt(p*(1-p)/float(n))

    df_possible = pd.read_csv('{0}/possible_counts_by_context_methyl.txt'.format(output_dir), 
                        sep='\t', index_col=['context','ref','alt','methylation_level']).rename(
        columns={'variant_count':'possible'})
    df_observed = pd.read_csv('{0}/observed_counts_by_context_methyl.txt'.format(output_dir), 
                        sep='\t', index_col=['context','ref','alt','methylation_level']).rename(
        columns={'variant_count':'observed'})
    df_po = df_possible.join(df_observed)
    df_po['proportion_observed'] = df_po['observed']/df_po['possible']

    df_mu = pd.read_csv('{0}/mu_by_context_methyl_downsampled_1000.txt'.format(output_dir),
                       sep='\t').rename(columns={'methyl_level':'methylation_level'})
    df_mu = df_mu.set_index(['context','ref','alt','methylation_level'])
    df_mu['sem'] = df_mu.apply(lambda row : sem(row['observed_1kg'], row['possible_variants']), axis = 1)

    df_po = df_po.join(df_mu[['mu','sem']])
    A, B = np.polyfit(df_po['mu'], np.log(1-df_po['proportion_observed']), 1,  w=1/df_po['sem'])
    df_po['fitted_po'] = 1-(np.exp(B))*(np.exp(A*df_po['mu']))
    print (r2_score(df_po['proportion_observed'], df_po['fitted_po'], sample_weight=1/df_po['sem']), (A,B))
    # 0.9987386340206728 (-18849750.846093018, -7.32454948743116e-05)

    # ouput mutation rate table
    df_po.drop(columns=['sem']).reset_index().to_csv(
        '{0}/mutation_rate_by_context_methyl.txt'.format(output_dir),
        sep='\t', quoting=csv.QUOTE_NONE, header=True, index=False)

    os.system('gsutil cp {0}/mutation_rate_by_context_methyl.txt {1}'.format(output_dir,output_bucket))


    ### Compute expected number of variants across the genome (per 1kb) based on the context-specific mutation rates

    # annotate by 1kb element
    context_ht = annotate_genome_element(context_ht,'{0}/misc/hg38.chrom.1kb.bed'.format(input_bucket))
    genome_ht = annotate_genome_element(genome_ht,'{0}/misc/hg38.chrom.1kb.bed'.format(input_bucket))    

    # count possible variants by context and methyl, per 1kb
    grouping = hl.struct(
        context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, methylation_level=context_ht.methyl_level, 
        element_id = context_ht.element_id,)
    output = {'variant_count': hl.agg.count()}
    possible_ht = context_ht.group_by(**grouping).aggregate(**output)
    possible_ht = possible_ht.persist()
    possible_ht.write('{0}/possible_counts_by_context_methyl_genome_1kb.ht'.format(output_bucket))
    # hl.read_table('{0}/possible_counts_by_context_methyl_genome_1kb.ht'.format(output_bucket)).export(
    #     '{0}/possible_counts_by_context_methyl_genome_1kb.txt'.format(output_bucket))

    # count observed variants per 1kb (save for computing contraint z scores)
    grouping = hl.struct(element_id = genome_ht.element_id)
    output = {'variant_count': hl.agg.count()}
    observed_ht = genome_ht.group_by(**grouping).aggregate(**output)
    observed_ht.write('{0}/observed_counts_genome_1kb.ht'.format(output_bucket))
    hl.read_table('{0}/observed_counts_genome_1kb.ht'.format(output_bucket)).export(
        '{0}/observed_counts_genome_1kb.txt'.format(output_bucket))

    # compute expected number of variants per 1kb
    po_ht = hl.import_table('{0}/mutation_rate_by_context_methyl.txt'.format(output_bucket),
                            delimiter='\t', 
                            types={'methylation_level':hl.tint, 
                                   'possible':hl.tint64, 
                                   'observed':hl.tint64,
                                   'proportion_observed':hl.tfloat,
                                   'mu':hl.tfloat,
                                   'fitted_po':hl.tfloat}).key_by('context','ref','alt','methylation_level')
    possible_ht = possible_ht.key_by('context','ref','alt','methylation_level')
    possible_ht = possible_ht.annotate(expected=possible_ht.variant_count*po_ht[possible_ht.key].fitted_po)
    possible_ht = possible_ht.key_by('element_id')

    expected_ht = possible_ht.group_by(possible_ht.element_id).aggregate(
        possible=hl.agg.sum(possible_ht.variant_count),
        expected=hl.agg.sum(possible_ht.expected))
    expected_ht.write('{0}/expected_counts_by_context_methyl_genome_1kb.ht'.format(output_bucket))
    expected_ht = hl.read_table('{0}/expected_counts_by_context_methyl_genome_1kb.ht'.format(output_bucket))
    expected_ht.annotate(expected = hl.format('%.8f', expected_ht.expected)).export(
        '{0}/expected_counts_by_context_methyl_genome_1kb.txt'.format(output_bucket))


    ### Train local sequence context and regional genomic features to predict DNM, per 1M

    # compute context-specific mutation rate for DNMs
    context_ht = hl.read_table('{0}/context_prefiltered.ht'.format(output_bucket))
    context_ht = filter_black_regions(context_ht)  

    dnm_ht = hl.import_table('{0}/misc/DNM_decode_psychencode.flip2hl.txt'.format(input_bucket),
                             delimiter='\t', 
                             types={'locus': hl.tlocus('GRCh38'), 'ref': hl.tstr, 'alt': hl.tstr}).key_by('locus')
    dnm_ht = dnm_ht.annotate(alleles = [dnm_ht.ref, dnm_ht.alt]).key_by('locus','alleles').select()
    dnm_ht = dnm_ht.repartition(2000)
    dnm_ht = context_ht.semi_join(dnm_ht)

    grouping = hl.struct(context=dnm_ht.context, ref=dnm_ht.ref, alt=dnm_ht.alt, methylation_level=dnm_ht.methyl_level)
    output = {'variant_count': hl.agg.count()}
    observed_ht = dnm_ht.group_by(**grouping).aggregate(**output)
    # observed_ht.write('{0}/observed_counts_by_context_methyl_dnm.ht'.format(output_bucket))

    grouping = hl.struct(context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, methylation_level=context_ht.methyl_level)
    output = {'variant_count': hl.agg.count()}
    possible_ht = context_ht.group_by(**grouping).aggregate(**output)
    # possible_ht.write('{0}/possible_counts_by_context_methyl_dnm.ht'.format(output_bucket))

    po_ht = possible_ht.annotate(
        observed = observed_ht[possible_ht.key].variant_count,
        proportion_observed = observed_ht[possible_ht.key].variant_count/possible_ht.variant_count)
    po_ht.rename({'variant_count' : 'possible'}).write(
        '{0}/proportion_observed_by_context_methyl_dnm.ht'.format(output_bucket))

    # count possible variants by context and methyl, per 1M
    context_ht = annotate_genome_element(context_ht,'{0}/misc/hg38.chrom.1M.bed'.format(input_bucket))
    grouping = hl.struct(
        context=context_ht.context, ref=context_ht.ref, alt=context_ht.alt, methylation_level=context_ht.methyl_level, 
        element_id = context_ht.element_id,)
    output = {'variant_count': hl.agg.count()}
    possible_ht = context_ht.group_by(**grouping).aggregate(**output)
    possible_ht.write('{0}/possible_counts_by_context_methyl_dnm_1M.ht'.format(output_bucket))

    # count observed DNMs per 1M (save for training)
    dnm_ht = annotate_genome_element(dnm_ht,'{0}/misc/hg38.chrom.1M.bed'.format(input_bucket))
    grouping = hl.struct(element_id = dnm_ht.element_id)
    output = {'variant_count': hl.agg.count()}
    observed_ht = dnm_ht.group_by(**grouping).aggregate(**output)
    observed_ht.write('{0}/observed_counts_dnm_1M.ht'.format(output_bucket))
    hl.read_table('{0}/observed_counts_dnm_1M.ht'.format(output_bucket)).export(
        '{0}/observed_counts_dnm_1M.txt'.format(output_bucket))

    # compute expected DNMs per 1M
    possible_ht = hl.read_table('{0}/possible_counts_by_context_methyl_dnm_1M.ht'.format(output_bucket)).key_by(
        'context','ref','alt','methylation_level')
    po_ht = hl.read_table('{0}/proportion_observed_by_context_methyl_dnm.ht'.format(output_bucket)).key_by(
      'context','ref','alt','methylation_level')
    possible_ht = possible_ht.annotate(expected=possible_ht.variant_count*po_ht[possible_ht.key].proportion_observed)
    possible_ht = possible_ht.key_by('element_id')

    expected_ht = possible_ht.group_by(possible_ht.element_id).aggregate(
        possible=hl.agg.sum(possible_ht.variant_count),
        expected=hl.agg.sum(possible_ht.expected))
    expected_ht.write('{0}/expected_counts_by_context_methyl_dnm_1M.ht'.format(output_bucket))
    expected_ht = hl.read_table('{0}/expected_counts_by_context_methyl_dnm_1M.ht'.format(output_bucket))
    expected_ht.annotate(expected = hl.format('%.8f', expected_ht.expected)).export(
        '{0}/expected_counts_by_context_methyl_dnm_1M.txt'.format(output_bucket))

    os.system('gsutil cp {0}/expected_counts_by_context_methyl_dnm_1M.txt {1}'.format(output_bucket,output_dir))
    os.system('gsutil cp {0}/observed_counts_dnm_1M.txt {1}'.format(output_bucket,output_dir))

    # combine expected DNMs with regional genomic features, and train on the observed DNMs
    os.system('gsutil cp {0}/misc/genomic_features17_1M.txt {1}/tmp/'.format(input_bucket,output_dir))

    df_expected = pd.read_csv('{0}/expected_counts_by_context_methyl_dnm_1M.txt'.format(output_dir), sep='\t', index_col='element_id')
    df_ft = pd.read_csv('{0}/tmp/genomic_features17_1M.txt'.format(output_dir), sep='\t', index_col='element_id')
    df_x_order = ['dist2telo','dist2cent','GC_content','expected','RT_BG02','LCR','SINE','LINE','recomb_male','recomb_female',
                  'met_sperm','met_oocyte','met_preimplantation','met_pgc','Nucleosome','cDNM_maternal_05M','cDNM_paternal_05M','CpG_island']

    df_x = df_expected.join(df_ft, how='inner')[df_x_order]
    df_y = pd.read_csv('{0}/observed_counts_dnm_1M.txt'.format(output_dir), sep='\t', index_col='element_id')

    # random forest regression: observed ~ expected (by sequence context)+ genomic features
    df_xy = df_y.join([df_x], how = 'inner')
    df_x = df_x[df_x.index.isin(df_xy.index)].sort_index()
    df_y = df_y[df_y.index.isin(df_xy.index)].sort_index()

    kf = KFold(n_splits=10, shuffle=True, random_state=0)
    test_index = df_y.iloc[list(kf.split(df_y))[-1][1]].index.tolist()

    x_train = df_x[~df_x.index.isin(test_index)]
    x_test = df_x[df_x.index.isin(test_index)]
    y_train = df_y[~df_y.index.isin(test_index)]
    y_test = df_y[df_y.index.isin(test_index)]

    rf = RandomForestRegressor(n_estimators=300, n_jobs=-1, max_features=3/4,verbose=0, random_state=0)
    rf.fit(x_train, y_train)
    r2 = rf.score(x_test, y_test)
    y_pred = rf.predict(x_test)

    ### Apply RF model to the gnomAD dataset, per 1kb

    os.system('gsutil cp {0}/expected_counts_by_context_methyl_genome_1kb.txt {1}'.format(output_bucket,output_dir))
    os.system('gsutil cp {0}/observed_counts_genome_1kb.txt {1}'.format(output_bucket,output_dir))
    os.system('gsutil cp {0}/misc/genomic_features17_1kb.txt {1}/tmp/'.format(input_bucket,output_dir))
    os.system('gsutil cp {0}/misc/RF_f18_dnm_1M.pkl {1}/tmp/'.format(input_bucket,output_dir))

    # load features
    df_expected = pd.read_csv('{0}/expected_counts_by_context_methyl_genome_1kb.txt'.format(output_dir), sep='\t', index_col='element_id')
    df_ft = pd.read_csv('{0}/tmp/genomic_features17_1kb.txt'.format(output_dir), sep='\t', index_col='element_id')
    df_x_order = ['dist2telo','dist2cent','GC_content','expected','RT_BG02','LCR','SINE','LINE','recomb_male','recomb_female',
                  'met_sperm','met_oocyte','met_preimplantation','met_pgc','Nucleosome','cDNM_maternal_05M','cDNM_paternal_05M','CpG_island']
    df_x = df_expected.join(df_ft, how='inner')[df_x_order]

    # predict expected number of variants from sequence context + genomic features 
    rf = pickle.load(open('{0}/tmp/RF_f18_dnm_1M.pkl'.format(output_dir), 'rb'))
    y_pred = rf.predict(df_x)
    df_x['predicted'] = y_pred

    # compute z scores from observed and expected counts
    df_observed = pd.read_csv('{0}/observed_counts_genome_1kb.txt'.format(output_dir), sep='\t', index_col='element_id')
    df_z = df_expected[['possible']].join(
        df_x[['predicted']].rename(columns={'predicted':'expected'})).join(
        df_observed.rename(columns={'variant_count':'observed'}))
    df_z['observed'] = df_z['observed'].fillna(value=0)

    df_z['oe'] = df_z['observed']/df_z['expected']
    df_z['chisq'] = (df_z['observed']-df_z['expected'])**2 / df_z['expected']
    df_z['z'] = np.where(df_z['oe'] >= 1., (-1)*np.sqrt(df_z['chisq']), np.sqrt(df_z['chisq']))
    df_z = df_z[df_z['z'].between(-10,10)].drop(columns=['chisq']).dropna()
    df_z.reset_index().to_csv(
        '{0}/constraint_z_genome_1kb_unfiltered.txt'.format(output_dir), sep='\t', quoting=csv.QUOTE_NONE, header=True, index=False)

    # qc elements by pass_filters and coverage
    os.system('gsutil cp {0}/misc/genome_1kb_gnomad_v31_pass.txt {1}/tmp/'.format(input_bucket,output_dir))
    os.system('gsutil cp {0}/misc/genome_1kb_gnomad_v31_coverage.txt {1}/tmp/'.format(input_bucket,output_dir))
    os.system('gsutil cp {0}/misc/genome_1kb_coding_exons.txt {1}/tmp/'.format(input_bucket,output_dir))

    d_pass = dict([line.strip().split('\t')[0],float(line.strip().split('\t')[1])] for line in 
                  open('{0}/tmp/genome_1kb_gnomad_v31_pass.txt'.format(output_dir)).readlines())
    d_coverage = dict([line.strip().split('\t')[0],float(line.strip().split('\t')[1])] for line in 
                      open('{0}/tmp/genome_1kb_gnomad_v31_coverage.txt'.format(output_dir)).readlines())
    df_z['pass'] = df_z.index.map(d_pass)
    df_z['coverage'] = df_z.index.map(d_coverage)
    df_z['pass_qc'] = (df_z['pass']>=0.8) & (df_z['coverage'].between(25,35)) & (df_z['possible']>=1000)

    # annotate proportion of coding bases
    d_coding = dict([line.strip().split('\t')[0],float(line.strip().split('\t')[1])] for line in 
                      open('{0}/tmp/genome_1kb_coding_exons.txt'.format(output_dir)).readlines())
    df_z['coding_prop'] = df_z.index.map(d_coding)
    df_z.reset_index().drop(columns=['pass','coverage']).to_csv(
        '{0}/constraint_z_genome_1kb.txt'.format(output_dir), sep='\t', quoting=csv.QUOTE_NONE, header=True, index=False)

    os.system('gsutil cp {0}/constraint_z_genome_1kb.txt {1}'.format(output_dir, output_bucket))

    # delete tmp files
    os.system('mv generic.py {0}/tmp/'.format(output_dir))
    os.system('mv constraint_basics.py {0}/tmp/'.format(output_dir))
    os.system('mv nc_constraint_utils.py {0}/tmp/'.format(output_dir))
    os.system('mv {0}/*1M.txt {0}/tmp/'.format(output_dir))
    # only keep non-redundant files in local dir
    os.system('mv {0}/possible*.txt {0}/tmp/'.format(output_dir))
    os.system('mv {0}/observed*.txt {0}/tmp/'.format(output_dir))
    os.system('rm -rf {0}/tmp/'.format(output_dir))
    os.system('rm {0}/hail*log'.format(output_dir))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-output_bucket', help='output gs://bucket_name/', required=True)
    parser.add_argument('-output_dir', help='output /path/to/local/dir_name', required=True)
    # parser.add_argument('-skip_hail', help='skip processing steps by Hail')
    args = parser.parse_args()
    main(args)




