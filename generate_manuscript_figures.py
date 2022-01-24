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
    plt_hist_freq_z("Figure_1a.pdf")
    print ("Plotting Fig.1b ...")
    plt_aps_vs_z("Figure_1b.pdf")     

    print ("Plotting Fig.2a ...")
    plt_enrichment_re("Figure_2a.pdf")
    print ("Plotting Fig.2b ...")
    plt_enrichment_9p21("Figure_2b.pdf")
    print ("Plotting Fig.2c ...")
    plt_9p21_locus("Figure_2c.pdf")

    print ("Plotting Fig.3a ...")
    plt_enrichment_gwas("Figure_3a.pdf")
    print ("Plotting Fig.3b ...")
    plt_enrichment_gwas_ukb("Figure_3b.pdf")

    print ("Plotting Fig.4a ...")
    plt_dd_cnv_z("Figure_4a.pdf")
    print ("Plotting Fig.4b ...")
    plt_dd_cnv_logit("Figure_4b.pdf")
    print ("Plotting Fig.4c ...")
    plt_cnv_ihh("Figure_4c.pdf")

    print ("Plotting Fig.5a ...")
    plt_prop_roadmaplinks("Figure_5a.pdf")
    print ("Plotting Fig.5b ...")
    plt_enh_geneset_z("Figure_5b.pdf")
    print ("Plotting Fig.5c ...")
    plt_geneset_loeuf0_enhz1("Figure_5c.pdf")
    print ("Plotting Fig.5d ...")
    plt_phastcons_loeuf0_enhz1("Figure_5d.pdf")
    print ("Plotting Fig.5e ...")
    plt_enhz_tissue_expr_corr("Figure_5e.pdf")

elif fig_x == "1":
    print ("Plotting Fig.1a ...")
    plt_hist_freq_z("Figure_1a.pdf")
    print ("Plotting Fig.1b ...")
    plt_aps_vs_z("Figure_1b.pdf")     
elif fig_x == "2":
    print ("Plotting Fig.2a ...")
    plt_enrichment_re("Figure_2a.pdf")
    print ("Plotting Fig.2b ...")
    plt_enrichment_9p21("Figure_2b.pdf")
    print ("Plotting Fig.2c ...")
    plt_9p21_locus("Figure_2c.pdf")
elif fig_x == "3":
    print ("Plotting Fig.3a ...")
    plt_enrichment_gwas("Figure_3a.pdf")
    print ("Plotting Fig.3b ...")
    plt_enrichment_gwas_ukb("Figure_3b.pdf")
elif fig_x == "4":
    print ("Plotting Fig.4a ...")
    plt_dd_cnv_z("Figure_4a.pdf")
    print ("Plotting Fig.4b ...")
    plt_dd_cnv_logit("Figure_4b.pdf")
    print ("Plotting Fig.4c ...")
    plt_cnv_ihh("Figure_4c.pdf")
elif fig_x == "5":
    print ("Plotting Fig.5a ...")
    plt_prop_roadmaplinks("Figure_5a.pdf")
    print ("Plotting Fig.5b ...")
    plt_enh_geneset_z("Figure_5b.pdf")
    print ("Plotting Fig.5c ...")
    plt_geneset_loeuf0_enhz1("Figure_5c.pdf")
    print ("Plotting Fig.5d ...")
    plt_phastcons_loeuf0_enhz1("Figure_5d.pdf")
    print ("Plotting Fig.5e ...")
    plt_enhz_tissue_expr_corr("Figure_5e.pdf")


