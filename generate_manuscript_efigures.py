#!/usr/bin/env python
# coding: utf-8

import argparse
import os
import sys

cwd = os.getcwd()
sys.path.append(cwd)
from efig_utils import *

parser = argparse.ArgumentParser()
parser.add_argument('-efig', help='figure number: 1, 2, 3, 4, 5, 6, or all', required=True)
args = parser.parse_args()

fig_x = args.efig

if fig_x.lower() == "all":

    print ("Plotting eFig.1a ...")
    plt_po_mu("eFigure_1a.pdf")
    print ("Plotting eFig.1b ...")
    plt_po_mu_chrx("eFigure_1b.pdf")
    print ("Plotting eFig.1c ...")
    plt_pred_dnms_ft("eFigure_1c.pdf")
    print ("Plotting eFig.1d ...")
    plt_scatter_oe_z("eFigure_1d.pdf")
    print ("Plotting eFig.1e ...")
    plt_scatter_oe_z_chrx("eFigure_1e.pdf")

    print ("Plotting eFig.2a ...")
    plt_hist_freq_z_exonic("eFigure_2a.pdf")
    print ("Plotting eFig.2b ...")
    plt_prop_z4_cd("eFigure_2b.pdf")
    print ("Plotting eFig.2c ...")
    plt_prop_z4_cd_titration("eFigure_2c.pdf")
    print ("Plotting eFig.2d ...")
    plt_z_ptcl_nc_cd("eFigure_2d.pdf")

    print ("Plotting eFig.3a ...")
    plt_enrichment_gwas_exl_ccre("eFigure_3a.pdf")
    print ("Plotting eFig.3b ...")
    plt_enrichment_gwas_ukb_exl_ccre("eFigure_3b.pdf")
    print ("Plotting eFig.3c ...")
    plt_gwas_ccre_z("eFigure_3c.pdf")

    print ("Plotting eFig.4a ...")
    plt_enrichment_re_10kb("eFigure_4a.pdf")
    print ("Plotting eFig.4b ...")
    plt_enrichment_gwas_10kb("eFigure_4b.pdf")
    print ("Plotting eFig.4c ...")
    plt_enrichment_gwas_ukb_10kb("eFigure_4c.pdf")

    print ("Plotting eFig.5a ...")
    plt_mave_enh_scores("eFigure_5a.pdf")
    print ("Plotting eFig.5b ...")
    plt_mave_enh_scores_corr("eFigure_5b.pdf")

    print ("Plotting eFig.6a ...")
    plt_hist_freq_z_chrx("eFigure_6a.pdf")
    print ("Plotting eFig.6b ...")
    plt_enrichment_re_chrx("eFigure_6b.pdf")


elif fig_x == "1":
    print ("Plotting eFig.1a ...")
    plt_po_mu("eFigure_1a.pdf")
    print ("Plotting eFig.1b ...")
    plt_po_mu_chrx("eFigure_1b.pdf")
    print ("Plotting eFig.1c ...")
    plt_pred_dnms_ft("eFigure_1c.pdf")
    print ("Plotting eFig.1d ...")
    plt_scatter_oe_z("eFigure_1d.pdf")
    print ("Plotting eFig.1e ...")
    plt_scatter_oe_z_chrx("eFigure_1e.pdf")
elif fig_x == "2":
    print ("Plotting eFig.2a ...")
    plt_hist_freq_z_exonic("eFigure_2a.pdf")
    print ("Plotting eFig.2b ...")
    plt_prop_z4_cd("eFigure_2b.pdf")
    print ("Plotting eFig.2c ...")
    plt_prop_z4_cd_titration("eFigure_2c.pdf")
    print ("Plotting eFig.2d ...")
    plt_z_ptcl_nc_cd("eFigure_2d.pdf")
elif fig_x == "3":
    print ("Plotting eFig.3a ...")
    plt_enrichment_gwas_exl_ccre("eFigure_3a.pdf")
    print ("Plotting eFig.3b ...")
    plt_enrichment_gwas_ukb_exl_ccre("eFigure_3b.pdf")
    print ("Plotting eFig.3c ...")
    plt_gwas_ccre_z("eFigure_3c.pdf")
elif fig_x == "4":
    print ("Plotting eFig.4a ...")
    plt_enrichment_re_10kb("eFigure_4a.pdf")
    print ("Plotting eFig.4b ...")
    plt_enrichment_gwas_10kb("eFigure_4b.pdf")
    print ("Plotting eFig.4c ...")
    plt_enrichment_gwas_ukb_10kb("eFigure_4c.pdf")
elif fig_x == "5":
    print ("Plotting eFig.5a ...")
    plt_mave_enh_scores("eFigure_5a.pdf")
    print ("Plotting eFig.5b ...")
    plt_mave_enh_scores_corr("eFigure_5b.pdf")
elif fig_x == "6":
    print ("Plotting eFig.6a ...")
    plt_hist_freq_z_chrx("eFigure_6a.pdf")
    print ("Plotting eFig.6b ...")
    plt_enrichment_re_chrx("eFigure_6b.pdf")





    