'''
Author: Taylor Reiter
Affiliation: UC Davis Lab for Data Intensive Biology
Aim: A Snakemake workflow to analyze differences between bins & NBHDs (donuts & crumbs) 
Date: May 4 2018
Run: snakemake --use-conda
Build graphic: snakemake --dag | dot -Tpdf > dag.pdf
'''

# Add installs:
# conda install -c bioconda bioconductor-biostrings

include: "crumb-bin.snakefile"
#include: "crumb-gold.snakefile"
include: "crumb-bin-kegg.snakefile"
#include: "crumb-gold-kegg.snakefile"

rule all:
    input:
        'outputs/hu-crumbs_bin/summary_of_inputs.tsv',
        dynamic('outputs/hu-crumbs_gold/unitigs/prokka-all/{crumb_gold}.faa'),
        #dynamic('outputs/hu-crumbs_gold/subtracts/prokka-all/{crumb_gold}.faa'),
        dynamic('outputs/hu-crumbs_gold/assembly/prokka/{crumb_gold}.faa'),
        # KEGG targets        
        #dynamic('outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/{crumb_bin}.pdf'),
        #dynamic('outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/{crumb_bin}.pdf'),
        #'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        #'outputs/hu-crumbs_bin/unitigs/an-et-al/an_uni.tsv',
        #'outputs/hu-crumbs_bin/unitigs/an-et-al/an_uni.pdf',
        #'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        #'outputs/hu-crumbs_bin/assembly/an-et-al/an_ass.tsv',
        #'outputs/hu-crumbs_bin/assembly/an-et-al/an_ass.pdf',
        #'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        #'outputs/hu-crumbs_bin/subtracts/an-et-al/an_sub.tsv',
        #'outputs/hu-crumbs_bin/subtracts/an-et-al/an_sub.pdf'

        # Gold targets
        #'outputs/hu-crumbs_gold/summary_of_inputs.tsv',
        #dynamic('outputs/hu-crumbs_gold/unitigs/prokka-all/{crumb_gold}.faa'),
        #dynamic('outputs/hu-crumbs_gold/subtracts/prokka-all/{crumb_gold}.faa'),
        #dynamic('outputs/hu-crumbs_gold/assembly/prokka/{crumb_gold}.faa'),
        #dynamic('outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/{crumb_gold}.pdf'),
        #dynamic('outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/{crumb_gold}.pdf'),
        #'outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        #'outputs/hu-crumbs_gold/unitigs/an-et-al/an_uni.tsv',
        #'outputs/hu-crumbs_gold/unitigs/an-et-al/an_uni.pdf',
        #'outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        #'outputs/hu-crumbs_gold/assembly/an-et-al/an_ass.tsv',
        #'outputs/hu-crumbs_gold/assembly/an-et-al/an_ass.pdf',
        #'outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        #'outputs/hu-crumbs_gold/subtracts/an-et-al/an_sub.tsv',
        #'outputs/hu-crumbs_gold/subtracts/an-et-al/an_sub.pdf'
