from snakemake.utils import R

# 1. 0.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa
# 2. 0.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa
# 3. 0.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa

# 1 is the unitigs from the crumb_bin
# 2 is #1 - #3
# 3 is the assembly of the crumb_bin reads using megahit

rule download_crumbs_bin:
    output: 'inputs/hu-crumbs_bin/hu-crumbs_bin.tar.gz'
    shell:'''
    curl -L -o {output} https://osf.io/ng4rw/download
    '''

rule unpack_crumbs_bin:
    output: 
        dynamic('inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa'),
        dynamic('inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa'),
        dynamic('inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa')
    input: 'inputs/hu-crumbs_bin/hu-crumbs_bin.tar.gz'
    params: output_folder = 'inputs/hu-crumbs_bin/'
    shell:'''
    tar xf {input} --directory {params.output_folder}
    '''

# SUMMARY STATS ----------------------------------------------------------

rule index_hu_unitigs:
    output: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.fai'
    input: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
    
rule index_hu_assembly:
    output: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa.fai'
    input: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
   
rule index_hu_subtract:
    output: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa.fai'
    input: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
    
rule summarize_hu:
    output: 'outputs/hu-crumbs_bin/summary_of_inputs.tsv'
    input: 
        dynamic('inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.fai'),
        dynamic('inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa.fai'),
        dynamic('inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa.fai')
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/skim_input.R
    '''

# UNITIGS ################################################################   
# megahit & annotate unitigs ---------------------------------------------
    
rule assemble_crumb_bin_unitigs:
    output: 'outputs/hu-crumbs_bin/unitigs/megahit/{crumb_bin}.contigs.fa'
    input: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/unitigs/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.crumb_bin} --out-prefix {wildcards.crumb_bin} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.crumb_bin}/{wildcards.crumb_bin}.contigs.fa {params.output_folder}/{wildcards.crumb_bin}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.crumb_bin}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_crumb_bin_unitigs:
    output: 'outputs/hu-crumbs_bin/unitigs/megahit-prokka/{crumb_bin}.faa'
    input:  'outputs/hu-crumbs_bin/unitigs/megahit/{crumb_bin}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/unitigs/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}mhuni
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------

rule prokka_crumbs_bin_unitigs:
    output: 'outputs/hu-crumbs_bin/unitigs/unitig-prokka/{crumb_bin}.faa'
    input:  'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/unitigs/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}uni
    touch {output}
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_crumb_bin_unitigs:
    output: 'outputs/hu-crumbs_bin/unitigs/prokka-all/{crumb_bin}.faa'
    input:  
        mh = 'outputs/hu-crumbs_bin/unitigs/megahit-prokka/{crumb_bin}.faa',
        uni = 'outputs/hu-crumbs_bin/unitigs/unitig-prokka/{crumb_bin}.faa'
    shell:'''
    Rscript --vanilla scripts/merge_fasta.R {input.mh} {input.uni} {output}
    '''

# blast ------------------------------------------------------------------

rule blastp_crumb_bin_unitigs:
    output: 'outputs/hu-crumbs_bin/unitigs/blast/{crumb_bin}-blastp.tab'
    input: 
        db = 'inputs/blast_db/interesting-aa.faa.psq',
        query = 'outputs/hu-crumbs_bin/unitigs/prokka-all/{crumb_bin}.faa'
    conda: 'env.yml'   
    shell: '''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
    '''
    
# SUBTRACTION ################################################################

rule assemble_crumb_bin_subtracts:
    output: 'outputs/hu-crumbs_bin/subtracts/megahit/{crumb_bin}.contigs.fa'
    input: 'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/subtracts/megahit'
    shell:'''
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.crumb_bin} --out-prefix {wildcards.crumb_bin} 
    mv {wildcards.crumb_bin}/{wildcards.crumb_bin}.contigs.fa {params.output_folder}/{wildcards.crumb_bin}.contigs.fa
    rm -rf {wildcards.crumb_bin}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_crumbs_bin_subtracts:
    output: 'outputs/hu-crumbs_bin/subtracts/megahit-prokka/{crumb_bin}.faa'
    input:  'outputs/hu-crumbs_bin/subtracts/megahit/{crumb_bin}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/subtracts/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}mhsub
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------
    
rule prokka_unitig_crumbs_bin_subtracts:
    output: 'outputs/hu-crumbs_bin/subtracts/unitig-prokka/{crumb_bin}.faa'
    input:  'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa.sub.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/subtracts/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}sub
    touch {output} 
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_subtracts:
    output: 'outputs/hu-crumbs_bin/subtracts/prokka-all/{crumb_bin}.faa'
    input:  
        mh = 'outputs/hu-crumbs_bin/subtracts/megahit-prokka/{crumb_bin}.faa',
        uni = 'outputs/hu-crumbs_bin/subtracts/unitig-prokka/{crumb_bin}.faa'
    shell:'''
    Rscript --vanilla scripts/merge_fasta.R {input.mh} {input.uni} {output}
    '''

# ASSEMBLIES ################################################################

# annotate megahit assemblies and hu genomes with prokka
rule prokka_crumbs_bin_assemblies:
    output: 'outputs/hu-crumbs_bin/assembly/prokka/{crumb_bin}.faa'
    input:  'inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs_bin/assembly/prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}ass
    touch {output}
    '''
