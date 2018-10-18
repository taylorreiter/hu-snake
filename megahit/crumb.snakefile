from snakemake.utils import R

rule download_crumbs:
    output: 'inputs/hu-crumbs-bin.tar.gz'
    shell:'''
    curl -L -o {output} https://osf.io/yxfv7/download
    '''

rule unpack_crumbs:
    output: 
        dynamic('inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumbs.fa'),
        dynamic('inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa')
    input: 'inputs/hu-crumbs-bin.tar.gz'
    params: output_folder = 'inputs/hu-s1-crumbs-aug7/'
    shell:'''
    tar xf {input} --directory {params.output_folder}
    '''

# SUMMARY STATS ----------------------------------------------------------

rule index_hu_unitigs:
    output: 'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.fai'
    input: 'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
    
rule index_hu_assembly:
    output: 'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa.fai'
    input: 'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
   
rule summarize_hu:
    output: 'outputs/hu-crumbs-bin/summary_of_inputs.tsv'
    input: 
        dynamic('inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa.fai'),
        dynamic('inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumbs.fa.assembly.fa.fai')
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/skim_input.R
    '''

# UNITIGS ################################################################   
# megahit & annotate unitigs ---------------------------------------------
    
rule assemble_crumb_unitigs:
    output: 'outputs/hu-crumbs-bin/unitigs/megahit/{crumb_bin}.contigs.fa'
    input: 'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs-bin/unitigs/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.crumb_bin} --out-prefix {wildcards.crumb_bin} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.crumb_bin}/{wildcards.crumb_bin}.contigs.fa {params.output_folder}/{wildcards.crumb_bin}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.crumb_bin}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_crumb_unitigs:
    output: 'outputs/hu-crumbs-bin/unitigs/megahit-prokka/{crumb_bin}.faa'
    input:  'outputs/hu-crumbs-bin/unitigs/megahit/{crumb_bin}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs-bin/unitigs/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}mhuni
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------

rule prokka_crumb_unitigs:
    output: 'outputs/hu-crumbs-bin/unitigs/unitig-prokka/{crumb_bin}.faa'
    input:  'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumbs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs-bin/unitigs/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}uni
    touch {output}
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_crumb_unitigs:
    output: 'outputs/hu-crumbs-bin/unitigs/prokka-all/{crumb_bin}.faa'
    input:  
        mh = 'outputs/hu-crumbs-bin/unitigs/megahit-prokka/{crumb_bin}.faa',
        uni = 'outputs/hu-crumbs-bin/unitigs/unitig-prokka/{crumb_bin}.faa'
    conda: 'env.yml'
    shell:'''
    Rscript --vanilla scripts/merge_fasta.R {input.mh} {input.uni} {output}
    '''
    
# ASSEMBLIES ################################################################

# assemble reads
rule megahit_crumb_reads:
    output: 'outputs/hu-crumbs-bin/assembly/megahit/{crumb_bin}.contigs.fa'
    input: 'inputs/hu-s1-crumbs-aug7/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumbs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs-bin/assembly/megahit'
    shell:'''
    megahit -r {input} --min-contig-len 200 --out-dir {wildcards.crumb_bin} --out-prefix {wildcards.crumb_bin} 
    mv {wildcards.crumb_bin}/{wildcards.crumb_bin}.contigs.fa {params.output_folder}/{wildcards.crumb_bin}.contigs.fa
    rm -rf {wildcards.crumb_bin}
    '''

# annotate megahit assemblies and hu genomes with prokka
rule prokka_crumb_assemblies:
    output: 
        faa = 'outputs/hu-crumbs-bin/assembly/prokka/{crumb_bin}.faa',
        fna = 'outputs/hu-crumbs-bin/assembly/prokka/{crumb_bin}.ffn'
    input: 'outputs/hu-crumbs-bin/assembly/megahit/{crumb_bin}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-crumbs-bin/assembly/prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.crumb_bin} --metagenome --force --locustag {wildcards.crumb_bin}ass
    touch {output.faa}
    '''

# collapse prokka annotations for kegg upload
rule cat_faa_crumb_assemblies:
    output: 'outputs/hu-crumbs-bin/assembly/prokka/all.faa'
    input: dynamic('outputs/hu-crumbs-bin/assembly/prokka/{crumb_bin}.faa')
    shell:'''
    cat outputs/hu-crumbs-bin/assembly/prokka/*faa > {output}
    '''

rule cat_ffn_crumb_assemblies:
    output: 'outputs/hu-crumbs-bin/assembly/prokka/all.ffn'
    input: dynamic('outputs/hu-crumbs-bin/assembly/prokka/{crumb_bin}.ffn')
    shell:'''
    cat outputs/hu-crumbs-bin/assembly/prokka/*.ffn > {output}
    '''
