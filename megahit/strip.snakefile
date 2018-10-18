from snakemake.utils import R


# CRUMBS ######################################################################## 
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


# CRUMB ASSEMBLIES ################################################################

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

# BINS ##########################################################################

rule download_bins:
    output: 'inputs/hu-bins.tar.gz'
    shell:
        "curl -L -o {output} https://osf.io/ehgbv/download"

rule unpack_bins:
    output: dynamic('inputs/hu-bins/{bin}.fa')
    input: 'inputs/hu-bins.tar.gz' 
    params: output_folder = 'inputs/hu-bins/'
    shell:
        "tar xf {input} --directory {params.output_folder}"

rule prokka_bins:
    output: 'outputs/hu-bins/prokka/{bin}.faa'
    input:  'inputs/hu-bins/{bin}.fa'
    conda: 'env.yml'
    params: outdir = 'outputs/hu-bins/prokka'
    shell:"""
    prokka {input} --outdir {params.outdir} --prefix {wildcards.bin} --metagenome --force --locustag {wildcards.bin}
    touch {output}
    """

rule cat_bins_faa:
    output: 'outputs/hu-bins/prokka/all-bins.faa'
    input: dynamic('outputs/hu-bins/prokka/{bin}.faa')
    shell:
        "cat {input} > {output}"
