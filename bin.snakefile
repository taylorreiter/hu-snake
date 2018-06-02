# rule all:
#     input: 
#         expand('{outputDIR}/{sample}.faa', outputDIR = config['outputDirProkkaG'], sample = config['sampleLIST'])

# download either the megahit assembled genomes (done by titus), or the genomes assembled by Hu et al. 
# Here, I download the genomes assembled by hu et al as posted to osf by titus

configfile: "config.yml"

rule download_hu_bins:
    output: 'inputs/hu-bins/{bin}.fa'
    shell:'''
    curl -L -o hu-bins.tar.gz https://osf.io/ehgbv/download
    mv hu-bins.tar.gz inputs/hu-bins/hu-bins.tar.gz
    cd inputs/hu-bins
    tar xf hu-bins.tar.gz
    touch {output}
    '''

# Prokka -----------------------------------------------------------

rule prokka_bins:
    output: 'outputs/hu-bins/prokka/{bin}.faa'
    input:  'inputs/hu-bins/{bin}.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-bins/prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.bin} --metagenome --force
    touch {output}
    '''
    
# BUSCO -------------------------------------------------------------
# https://www.microbe.net/2017/12/13/why-genome-completeness-and-contamination-estimates-are-more-complicated-than-you-think/


    
rule run_busco_bac_bins:
    output: 'outputs/hu-bins/busco/run_{bin}_bac'
    input: 
        bin_in='inputs/hu-bins/{bin}.fa',
        busco_db='inputs/busco/bacteria_odb9/'
    conda:  "env.yml"
    shell:'''
	run_busco -i {input.bin_in} -o {wildcards.bin}_bac -l {input.busco_db} -m geno
    mv run_{wildcards.bin}_bac {output}
    '''

