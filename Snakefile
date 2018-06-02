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
include: "crumb-gold.snakefile"
include: "crumb-bin-kegg.snakefile"
include: "crumb-gold-kegg.snakefile"

rule all:
    input:
        #dynamic('outputs/hu-crumbs_bin/unitigs/blast/{crumb_bin}-blastp.tab'),
        #dynamic('outputs/hu-crumbs_bin/subtracts/blast/{crumb_bin}-blastp.tab'),
        #dynamic('outputs/hu-crumbs_bin/assembly/blast/{crumb_bin}-blastp.tab'),
        #dynamic('outputs/hu-crumbs_bin/unitigs/busco/run_{crumb_bin}_bac'),
        #dynamic('outputs/hu-crumbs_bin/assembly/busco/run_{crumb_bin}_bac'),
        'outputs/hu-crumbs_bin/summary_of_inputs.tsv',
        #dynamic('outputs/hu-crumbs_bin/assembly/blast/{crumb_bin}-tax-blastp.tab'),
        #dynamic('outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/{crumb_bin}.pdf'),
        #dynamic('outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/{crumb_bin}.pdf'),
        'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-each.pdf',
        'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        'outputs/hu-crumbs_bin/unitigs/an-et-al/an_uni.tsv',
        'outputs/hu-crumbs_bin/unitigs/an-et-al/an_uni.pdf',
        'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-each.pdf',
        'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        'outputs/hu-crumbs_bin/assembly/an-et-al/an_ass.tsv',
        'outputs/hu-crumbs_bin/assembly/an-et-al/an_ass.pdf',
        'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-each.pdf',
        'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        'outputs/hu-crumbs_bin/subtracts/an-et-al/an_sub.tsv',
        'outputs/hu-crumbs_bin/subtracts/an-et-al/an_sub.pdf'
        
