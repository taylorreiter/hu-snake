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
        # expand('outputs/hu-bins/prokka/{bin}.faa', bin = config['sampleBIN']),
        # expand('outputs/hu-bins/busco/run_{bin}_bac', bin =  config['sampleBIN']),
        # expand('outputs/hu-croissants/busco/run_{croissant}_bac', croissant = config['sampleCROISSANT']),
        dynamic('outputs/hu-croissants/unitigs/blast/{croissant}-blastp.tab'),
        dynamic('outputs/hu-croissants/subtracts/blast/{croissant}-blastp.tab'),
        dynamic('outputs/hu-croissants/assembly/blast/{croissant}-blastp.tab')
        