#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys

cwd = os.getcwd()
sys.path.append(cwd)
from efig_utils import *

parser = argparse.ArgumentParser()
parser.add_argument("-efig", help="figure number: 1, 2, 3, 4, 5, 6, or all", required=True)
args = parser.parse_args()

fig_x = args.efig

if fig_x.lower() == "all":

    print ("Plotting eFig.1a ...")
    plt_po_mu("eFigure_1a.pdf")
    print ("Plotting eFig.1b ...")
    plt_po_mu_chrx("eFigure_1b.pdf")
    print ("Plotting eFig.1c ...")
    plt_genomic_ft_sel("eFigure_1c.pdf")
    print ("Plotting eFig.1d ...")
    plt_scatter_oe_gnocchi("eFigure_1d.pdf")
    print ("Plotting eFig.1e ...")
    plt_scatter_oe_gnocchi_chrx("eFigure_1e.pdf")

    print ("Plotting eFig.2a ...")
    plt_prop_gnocchi4_cd("eFigure_2a.pdf")
    print ("Plotting eFig.2b ...")
    plt_hist_freq_gnocchi_exonic("eFigure_2b.pdf")
    print ("Plotting eFig.2c ...")
    plt_prop_gnocchi4_cd_titration("eFigure_2c.pdf")
    print ("Plotting eFig.2d ...")
    plt_gnochhi_ptcl_nc_cd("eFigure_2d.pdf")

    print ("Plotting eFig.3a ...")
    plt_gnocchi_ptcl_annot_cd("eFigure_3a.pdf")
    print ("Plotting eFig.3b ...")
    plt_gnocchi_ptcl_annot_nc("eFigure_3b.pdf")

    print ("Plotting eFig.4a ...")
    plt_enrichment_gwas_vs_ccre("eFigure_4a.pdf")
    print ("Plotting eFig.4b ...")
    plt_delta_pip("eFigure_4b.pdf")

    print ("Plotting eFig.5a ...")
    plt_comparison_roc("gwas_catalog","topmed_maf5",0,"eFigure_5a-1.pdf")
    plt_comparison_roc("gwas_fine-mapping","topmed_maf5",0,"eFigure_5a-2.pdf")
    plt_comparison_roc("gwas_fine-mapping_hc","topmed_maf5",0,"eFigure_5a-3.pdf")
    plt_comparison_roc("clinvar_plp_hgmd","topmed_mac1",0,"eFigure_5a-4.pdf")
    print ("Plotting eFig.5b ...")
    plt_comparison_auc_af("gwas_catalog","eFigure_5b-1.pdf")
    plt_comparison_auc_af("gwas_fine-mapping","eFigure_5b-2.pdf")
    plt_comparison_auc_af("gwas_fine-mapping_hc","eFigure_5b-3.pdf")
    plt_comparison_auc_af("clinvar_plp_hgmd","eFigure_5b-4.pdf")

    print ("Plotting eFig.6 ...")
    plt_comparison_roc_gnocchi("gwas_catalog","topmed_maf5",0,"eFigure_6-1.pdf")
    plt_comparison_roc_gnocchi("gwas_fine-mapping","topmed_maf5",0,"eFigure_6-2.pdf")
    plt_comparison_roc_gnocchi("gwas_fine-mapping_hc","topmed_maf5",0,"eFigure_6-3.pdf")
    plt_comparison_roc_gnocchi("clinvar_plp_hgmd","topmed_mac1",0,"eFigure_6-4.pdf")

    print ("Plotting eFig.7 ...")
    plt_score_corr(["Gnocchi","Orion","CDTS","gwRVIS","DR","phyloP","phastCons","GERP"],"eFigure_7.pdf")

    print ("Plotting eFig.8a ...")
    plt_power_depl(1000,"eFigure_8a.pdf")
    print ("Plotting eFig.8b ...")
    plt_power_depl(100,"eFigure_8b.pdf")
    print ("Plotting eFig.8c ...")
    plt_comparison_auc_ws("eFigure_8c.pdf")
    print ("Plotting eFig.8d ...")
    plt_comparison_auc_pop("eFigure_8d.pdf")
        

elif fig_x == "1":
    print ("Plotting eFig.1a ...")
    plt_po_mu("eFigure_1a.pdf")
    print ("Plotting eFig.1b ...")
    plt_po_mu_chrx("eFigure_1b.pdf")
    print ("Plotting eFig.1c ...")
    plt_genomic_ft_sel("eFigure_1c.pdf")
    print ("Plotting eFig.1d ...")
    plt_scatter_oe_gnocchi("eFigure_1d.pdf")
    print ("Plotting eFig.1e ...")
    plt_scatter_oe_gnocchi_chrx("eFigure_1e.pdf")
elif fig_x == "2":
    print ("Plotting eFig.2a ...")
    plt_prop_gnocchi4_cd("eFigure_2a.pdf")
    print ("Plotting eFig.2b ...")
    plt_hist_freq_gnocchi_exonic("eFigure_2b.pdf")
    print ("Plotting eFig.2c ...")
    plt_prop_gnocchi4_cd_titration("eFigure_2c.pdf")
    print ("Plotting eFig.2d ...")
    plt_gnochhi_ptcl_nc_cd("eFigure_2d.pdf")
elif fig_x == "3":
    print ("Plotting eFig.3a ...")
    plt_gnocchi_ptcl_annot_cd("eFigure_3a.pdf")
    print ("Plotting eFig.3b ...")
    plt_gnocchi_ptcl_annot_nc("eFigure_3b.pdf")
elif fig_x == "4":
    print ("Plotting eFig.4a ...")
    plt_enrichment_gwas_vs_ccre("eFigure_4a.pdf")
    print ("Plotting eFig.4b ...")
    plt_delta_pip("eFigure_4b.pdf")
elif fig_x == "5":
    print ("Plotting eFig.5a ...")
    plt_comparison_roc("gwas_catalog","topmed_maf5",0,"eFigure_5a-1.pdf")
    plt_comparison_roc("gwas_fine-mapping","topmed_maf5",0,"eFigure_5a-2.pdf")
    plt_comparison_roc("gwas_fine-mapping_hc","topmed_maf5",0,"eFigure_5a-3.pdf")
    plt_comparison_roc("clinvar_plp_hgmd","topmed_mac1",0,"eFigure_5a-4.pdf")
    print ("Plotting eFig.5b ...")
    plt_comparison_auc_af("gwas_catalog","eFigure_5b-1.pdf")
    plt_comparison_auc_af("gwas_fine-mapping","eFigure_5b-2.pdf")
    plt_comparison_auc_af("gwas_fine-mapping_hc","eFigure_5b-3.pdf")
    plt_comparison_auc_af("clinvar_plp_hgmd","eFigure_5b-4.pdf")
elif fig_x == "6":
    print ("Plotting eFig.6 ...")
    plt_comparison_roc_gnocchi("gwas_catalog","topmed_maf5",0,"eFigure_6-1.pdf")
    plt_comparison_roc_gnocchi("gwas_fine-mapping","topmed_maf5",0,"eFigure_6-2.pdf")
    plt_comparison_roc_gnocchi("gwas_fine-mapping_hc","topmed_maf5",0,"eFigure_6-3.pdf")
    plt_comparison_roc_gnocchi("clinvar_plp_hgmd","topmed_mac1",0,"eFigure_6-4.pdf")
elif fig_x == "7":
    print ("Plotting eFig.7 ...")
    plt_score_corr(["Gnocchi","Orion","CDTS","gwRVIS","DR","phyloP","phastCons","GERP"],"eFigure_7.pdf")

elif fig_x == "8":
    print ("Plotting eFig.8a ...")
    plt_power_depl(1000,"eFigure_8a.pdf")
    print ("Plotting eFig.8b ...")
    plt_power_depl(100,"eFigure_8b.pdf")
    print ("Plotting eFig.8c ...")
    plt_comparison_auc_ws("eFigure_8c.pdf")
    print ("Plotting eFig.8d ...")
    plt_comparison_auc_pop("eFigure_8d.pdf")





    
