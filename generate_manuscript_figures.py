#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys

cwd = os.getcwd()
sys.path.append(cwd)
from fig_utils import *

parser = argparse.ArgumentParser()
parser.add_argument('-fig', help='figure number: 1, 2, 3, 4, 5, or all', required=True)
args = parser.parse_args()

fig_x = args.fig

if fig_x.lower() == "all":

    print ("Plotting Fig.1a ...")
    plt_hist_freq_gnocchi("Figure_1a.pdf")
    print ("Plotting Fig.1b ...")
    plt_aps_vs_gnocchi("Figure_1b.pdf")     

    print ("Plotting Fig.2a ...")
    plt_enrichment_re("Figure_2a.pdf")
    print ("Plotting Fig.2b ...")
    plt_enrichment_gwas("Figure_2b.pdf")
    print ("Plotting Fig.2c ...")
    plt_enrichment_gwas_ukb("Figure_2c.pdf")
    print ("Plotting Fig.2d ...")
    plt_plg_gwas("Figure_2d.pdf")

    print ("Plotting Fig.3a ...")
    plt_comparison_roc("gwas_fine-mapping","topmed_maf5",0,"Figure_3a.pdf")
    print ("Plotting Fig.3b ...")
    plt_comparison_roc("clinvar_plp_hgmd","topmed_mac1",0,"Figure_3b.pdf")
    print ("Plotting Fig.3c ...")
    plt_dominance_scores("gwas","Figure_3c.pdf")
    print ("Plotting Fig.3d ...")
    plt_dominance_scores("clinvar_plp_hgmd","Figure_3d.pdf")

    print ("Plotting Fig.4a ...")
    plt_cnv_dd_gnocchi("Figure_4a.pdf")
    print ("Plotting Fig.4b ...")
    plt_dd_cnv_logit("Figure_4b.pdf")
    print ("Plotting Fig.4c ...")
    plt_cnv_ihh("Figure_4c.pdf")
    print ("Plotting Fig.4d ...")
    plt_cnv_recurrent("Figure_4d-1.pdf")
    plt_cnv_recurrent_gnocchi("Figure_4d-2.pdf")

    print ("Plotting Fig.5a ...")
    plt_prop_roadmaplinks("Figure_5a.pdf")
    print ("Plotting Fig.5b ...")
    plt_enh_geneset_gnocchi("Figure_5b.pdf")
    print ("Plotting Fig.5c ...")
    plt_enh_gnocchi_loeuf_roc("Figure_5c.pdf")
    print ("Plotting Fig.5d ...")
    plt_enh_gnocchi_tissue_expr_corr("Figure_5d.pdf")


elif fig_x == "1":
    print ("Plotting Fig.1a ...")
    plt_hist_freq_gnocchi("Figure_1a.pdf")
    print ("Plotting Fig.1b ...")
    plt_aps_vs_gnocchi("Figure_1b.pdf")      
elif fig_x == "2":
    print ("Plotting Fig.2a ...")
    plt_enrichment_re("Figure_2a.pdf")
    print ("Plotting Fig.2b ...")
    plt_enrichment_gwas("Figure_2b.pdf")
    print ("Plotting Fig.2c ...")
    plt_enrichment_gwas_ukb("Figure_2c.pdf")
    print ("Plotting Fig.2d ...")
    plt_plg_gwas("Figure_2d.pdf")
elif fig_x == "3":
    print ("Plotting Fig.3a ...")
    plt_comparison_roc("gwas_fine-mapping","topmed_maf5",0,"Figure_3a.pdf")
    print ("Plotting Fig.3b ...")
    plt_comparison_roc("clinvar_plp_hgmd","topmed_mac1",0,"Figure_3b.pdf")
    print ("Plotting Fig.3c ...")
    plt_dominance_scores("gwas","Figure_3c.pdf")
    print ("Plotting Fig.3d ...")
    plt_dominance_scores("clinvar_plp_hgmd","Figure_3d.pdf")
elif fig_x == "4":
    print ("Plotting Fig.4a ...")
    plt_cnv_dd_gnocchi("Figure_4a.pdf")
    print ("Plotting Fig.4b ...")
    plt_dd_cnv_logit("Figure_4b.pdf")
    print ("Plotting Fig.4c ...")
    plt_cnv_ihh("Figure_4c.pdf")
    print ("Plotting Fig.4d ...")
    plt_cnv_recurrent("Figure_4d-1.pdf")
    plt_cnv_recurrent_gnocchi("Figure_4d-2.pdf")
elif fig_x == "5":
    print ("Plotting Fig.5a ...")
    plt_prop_roadmaplinks("Figure_5a.pdf")
    print ("Plotting Fig.5b ...")
    plt_enh_geneset_gnocchi("Figure_5b.pdf")
    print ("Plotting Fig.5c ...")
    plt_enh_gnocchi_loeuf_roc("Figure_5c.pdf")
    print ("Plotting Fig.5d ...")
    plt_enh_gnocchi_tissue_expr_corr("Figure_5d.pdf")

