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
