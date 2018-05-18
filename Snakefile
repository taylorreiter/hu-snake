'''
Author: Taylor Reiter
Affiliation: UC Davis Lab for Data Intensive Biology
Aim: A Snakemake workflow to analyze differences between bins & NBHDs (donuts & croissants) 
Date: May 4 2018
Run: snakemake --use-conda
Build graphic: snakemake --dag | dot -Tpdf > dag.pdf
'''

# Add installs:
# conda install -c bioconda bioconductor-biostrings

include: "croissant.snakefile"
include: "genome.snakefile"

rule all:
    input:
        #dynamic('outputs/hu-croissants/unitigs/blast/{croissant}-blastp.tab'),
        #dynamic('outputs/hu-croissants/subtracts/blast/{croissant}-blastp.tab'),
        #dynamic('outputs/hu-croissants/assembly/blast/{croissant}-blastp.tab'),
        dynamic('outputs/hu-croissants/unitigs/busco/run_{croissant}_bac'),
        dynamic('outputs/hu-croissants/assembly/busco/run_{croissant}_bac')
        