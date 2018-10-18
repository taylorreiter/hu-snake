from snakemake.utils import R

# 1. hu-genome*_gold.fa.cdbg_ids.contigs.fa.gz.crumbs.fa
# 2. hu-genome*_gold.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa
# 3. hu-genome*_gold.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa
# 4. hu-genome*_gold.fa.cdbg_ids.reads.fa.gz.crumbs.fa

# 1 is the unitigs from the crumb_gold
# 2 is #1 - #3
# 3 is the assembly of the crumb_gold reads using megahit
# 4 is the reads that went into the assembly

rule download_crumbs_gold:
    output: 'inputs/hu-crumbs_gold/hu-crumbs_gold.tar.gz'
    shell:'''
    curl -L -o {output} https://osf.io/tpwhz/download
    '''

rule unpack_crumbs_gold:
    output: 
        dynamic('inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'),
        dynamic('inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa'),
        dynamic('inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa')
    input: 'inputs/hu-crumbs_gold/hu-crumbs_gold.tar.gz'
    params: output_folder = 'inputs/hu-crumbs_gold/'
    shell:'''
    tar xf {input} --directory {params.output_folder}
    '''
       
# SUMMARY STATS ----------------------------------------------------------

rule index_gold_unitigs:
    output: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.fai'
    input: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
    
rule index_gold_assembly:
    output: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa.fai'
    input: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
   
rule index_gold_subtract:
    output: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa.fai'
    input: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
    
rule summarize_gold:
    output: 'outputs/hu-crumbs_gold/summary_of_inputs.tsv'
    input: 
        dynamic('inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.fai'),
        dynamic('inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa.fai'),
        dynamic('inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa.fai')
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/skim_input_gold.R
    '''

# UNITIGS ################################################################   
# megahit & annotate unitigs ---------------------------------------------
    
rule assemble_crumb_gold_unitigs:
    output: 'outputs/hu-crumbs_gold/unitigs/megahit/{crumb_gold}.contigs.fa'
    input: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/unitigs/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.crumb_gold} --out-prefix {wildcards.crumb_gold} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.crumb_gold}/{wildcards.crumb_gold}.contigs.fa {params.output_folder}/{wildcards.crumb_gold}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.crumb_gold}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_crumb_gold_unitigs:
    output: 'outputs/hu-crumbs_gold/unitigs/megahit-prokka/{crumb_gold}.faa'
    input:  'outputs/hu-crumbs_gold/unitigs/megahit/{crumb_gold}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/unitigs/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_gold} --metagenome --force --locustag {wildcards.crumb_gold}mhuni
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------

rule prokka_crumbs_gold_unitigs:
    output: 'outputs/hu-crumbs_gold/unitigs/unitig-prokka/{crumb_gold}.faa'
    input:  'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/unitigs/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_gold} --metagenome --force --locustag {wildcards.crumb_gold}uni
    touch {output}
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_crumb_gold_unitigs:
    output: 'outputs/hu-crumbs_gold/unitigs/prokka-all/{crumb_gold}.faa'
    input:  
        mh = 'outputs/hu-crumbs_gold/unitigs/megahit-prokka/{crumb_gold}.faa',
        uni = 'outputs/hu-crumbs_gold/unitigs/unitig-prokka/{crumb_gold}.faa'
    shell:'''
    Rscript --vanilla scripts/merge_fasta.R {input.mh} {input.uni} {output}
    '''

# SUBTRACTION ################################################################

rule assemble_crumb_gold_subtracts:
    output: 'outputs/hu-crumbs_gold/subtracts/megahit/{crumb_gold}.contigs.fa'
    input: 'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/subtracts/megahit'
    shell:'''
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.crumb_gold} --out-prefix {wildcards.crumb_gold} 
    mv {wildcards.crumb_gold}/{wildcards.crumb_gold}.contigs.fa {params.output_folder}/{wildcards.crumb_gold}.contigs.fa
    rm -rf {wildcards.crumb_gold}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_crumbs_gold_subtracts:
    output: 'outputs/hu-crumbs_gold/subtracts/megahit-prokka/{crumb_gold}.faa'
    input:  'outputs/hu-crumbs_gold/subtracts/megahit/{crumb_gold}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/subtracts/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_gold} --metagenome --force --locustag {wildcards.crumb_gold}mhsub
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------
    
rule prokka_unitig_crumbs_gold_subtracts:
    output: 'outputs/hu-crumbs_gold/subtracts/unitig-prokka/{crumb_gold}.faa'
    input:  'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.sub.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/subtracts/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_gold} --metagenome --force --locustag {wildcards.crumb_gold}sub
    touch {output} 
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_gold_subtracts:
    output: 'outputs/hu-crumbs_gold/subtracts/prokka-all/{crumb_gold}.faa'
    input:  
        mh = 'outputs/hu-crumbs_gold/subtracts/megahit-prokka/{crumb_gold}.faa',
        uni = 'outputs/hu-crumbs_gold/subtracts/unitig-prokka/{crumb_gold}.faa'
    shell:'''
    Rscript --vanilla scripts/merge_fasta.R {input.mh} {input.uni} {output}
    '''

# ASSEMBLIES ################################################################

# annotate megahit assemblies and hu genomes with prokka
rule prokka_crumbs_gold_assemblies:
    output: 'outputs/hu-crumbs_gold/assembly/prokka/{crumb_gold}.faa'
    input:  'inputs/hu-crumbs_gold/{crumb_gold}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_gold/assembly/prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_gold} --metagenome --force --locustag {wildcards.crumb_gold}ass
    touch {output}
    '''
