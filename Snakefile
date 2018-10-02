'''
Author: Taylor Reiter
Affiliation: UC Davis Lab for Data Intensive Biology
Aim: A Snakemake workflow to analyze differences between bins & NBHDs (donuts & crumbs) 
Date: May 4 2018
Run: snakemake --use-conda
Build graphic: snakemake --dag | dot -Tpdf > dag.pdf
'''

include: "strip.snakefile"
include: "strip-kegg.snakefile"

rule all:
    input:
        'outputs/hu-crumbs-bin/assembly/prokka/all.faa',
        'outputs/hu-crumbs-bin/assembly/prokka/all.ffn',
        'outputs/hu-bins/prokka/all-bins.faa',
        dynamic('outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/{crumb_bin}.pdf'),
        'outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/kegg-each.pdf',
        'outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        'outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        'outputs/hu-crumbs-bin/assembly/an-et-al/an_ass.tsv',
        'outputs/hu-crumbs-bin/assembly/an-et-al/an_ass.pdf',
        'outputs/figures/fig5b.pdf',
        'outputs/figures/fig5b.png',
        'outputs/figures/fig5a.pdf',
        'outputs/figures/fig5a.png',
        'outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc.csv',
        'outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc_full.csv',
        'outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv',
        'outputs/hu-crumbs-bin/assembly/marker-genes/blast_aas_full.csv',
        'outputs/hu-crumbs-bin/assembly/marker-genes/novel_crumb_marks.csv',
        'outputs/hu-crumbs-bin/assembly/other/crumb_nitro_genes.csv'
