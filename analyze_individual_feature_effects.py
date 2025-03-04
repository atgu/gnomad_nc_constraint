import pandas as pd
import numpy as np
import scipy
from scipy import stats
import statsmodels.api as sm

# Download variant & annotation files
input_bucket = 'gs://gnomad-nc-constraint-v31-paper'
if not os.path.exists('{0}/genomic_features'.format(output_dir)): os.mkdir('{0}/genomic_features'.format(output_dir))
os.system('gsutil cp {0}/genomic_features/* {1}/genomic_features'.format(input_bucket,output_dir))

# de novo variants
df_dnm1 = pd.read_csv('{0}/genomic_features/DNM_decode_psychencode_site_context.mutation_rate.txt'.format(output_dir),sep='\t')  
df_ft_1 = pd.read_csv('{0}/genomic_features/genomic_features13_dnm1_flnk_1k-1M.txt'.format(output_dir),sep='\t').drop_duplicates()
df_dnm1 = df_dnm1.merge(df_ft_1.rename(columns={'element_id':'locus'}),how='left',on='locus')
# ‘non-mutated’ background
df_dnm0 = pd.read_csv('{0}/genomic_features/context_prefiltered_nonmutated-dnm_sites10xdnm.mutation_rate.txt'.format(output_dir),sep='\t')
df_dnm0 = df_dnm0[~df_dnm0['locus'].str.contains('chrX:')]
df_ft_0 = pd.read_csv('{0}/genomic_features/genomic_features13_dnm0_10x_flnk_1k-1M.txt'.format(output_dir),sep='\t').drop_duplicates()
df_dnm0 = df_dnm0.merge(df_ft_0.rename(columns={'element_id':'locus'}),how='left',on='locus')

ft_cols = ['dist2telo', 'dist2cent', 'LCR', 'SINE', 'LINE','GC_content', 
           'recomb_male', 'recomb_female', 'met_sperm','Nucleosome','CpG_island',
           'cDNM_maternal_05M', 'cDNM_paternal_05M']

dnm_po = pd.read_csv('{0}/mutation_rate_by_context_methyl.txt'.format(output_dir),sep='\t')
contexts = sorted(list(set(dnm_po['context'])))

of = open('{0}/genomic_features/dnm01_10x_ft_logit_regularized_coef_z_3mer_context_flnk_1k-1M.txt'.format(output_dir),'w')
of.write('\t'.join(['context','window','feature','coef','se','pval']) + '\n')
for context in contexts:
    df_1_ = df_dnm1[(df_dnm1['context']==context)]
    df_0_ = df_dnm0[(df_dnm0['context']==context)]
    for flnk in ['1k','10k','100k','1M']:
        for ft_ in ft_cols:
            ft = ft_+'_'+flnk

            df_1 = df_1_[['locus',ft]].dropna()
            df_1['group'] = 1
            df_0 = df_0_[['locus',ft]].dropna()
            df_0['group'] = 0
            df_01 = pd.concat([df_1,df_0])

            df_y = df_01[['group']]
            df_x = df_01[[ft]]
            df_x = df_x.apply(scipy.stats.zscore)

            try:
                logit = sm.Logit(df_y, sm.add_constant(df_x[[ft]], has_constant='add')).fit_regularized()
            except:
                coef, pval, lci, hci, se = [np.nan]*5
            else:
                coef = logit.params[[ft]][0]
                pval = logit.pvalues[[ft]][0]
                lci, hci = logit.conf_int().transpose()[ft].to_list()
                se = (hci-coef)/1.96
            of.write('\t'.join([context,flnk,ft_,str(coef),str(se),str(pval)]) + '\n')
of.close()


# Select significant features
df = pd.read_csv('{0}/genomic_features/dnm01_10x_ft_logit_regularized_coef_z_3mer_context_flnk_1k-1M.txt'.format(output_dir),sep="\t")
df["bonf"] = np.where(df["context"].isin(["ACG","CCG","GCG","TCG"]),0.05/4/8, 0.05/4/13)
df_sig = df[df["pval"]<=df["bonf"]]
df_sig = df_sig.drop_duplicates(subset=["context","feature"])
df_sig[["context","feature","window","coef","se","pval"]].to_csv(
    '{0}/genomic_features/dnm01_10x_ft_logit_regularized_coef_z_3mer_context_flnk_1k-1M.selected.txt'.format(output_dir),
           sep="\t", quoting=csv.QUOTE_NONE, header=True, index=False)
# this file corresponds to gs://gnomad-nc-constraint-v31-paper/misc/genomic_features13_sel.txt

