# download and assemble unitigs from nbhds.
rule download_donuts:
    output: 'inputs/hu-donuts/{sample}.fa.nbhd.fa.gz'
    shell:'''
    curl -L -o inputs/hu-donuts/hu-genome08.fa.nbhd.fa.gz https://osf.io/g2zvh/download 
    curl -L -o inputs/hu-donuts/hu-genome11.fa.nbhd.fa.gz https://osf.io/5k2zj/download
    curl -L -o inputs/hu-donuts/hu-genome19.fa.nbhd.fa.gz https://osf.io/29bk8/download
    ''' 
    
rule assemble_donuts:
    output: '{outputDIR}/{sample}.contigs.fa'
    input: 'inputs/hu-donuts/{sample}.fa.nbhd.fa.gz'
    conda: 'env.yml'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --out-dir {wildcards.sample} --out-prefix {wildcards.sample} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.sample}/{wildcards.sample}.contigs.fa {wildcards.outputDIR}/{wildcards.sample}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.sample}
    ''' 

rule prokka_donuts:
    output: '{outputDIR}/{sample}.faa'
    input:  '{outputDIR}/{sample}.contigs.fa'
    conda: 'env.yml'
    shell:'''
    prokka {input} --outdir {wildcards.outputDIR} --prefix {wildcards.sample} --metagenome --force
    touch {output}
    '''
    
# rule all:
#    input: 'outputs/busco/run_wild_olive_busco'
# 
# rule download_sylv_inputs_busco:
#     output: 
#         gz='inputs/sylvestris/Olea_europaea_all_scaffolds.fa.gz',
#         uncmp='inputs/sylvestris/Olea_europaea_all_scaffolds.fa'
#     shell:'''
# 	wget -O {output.gz} http://olivegenome.org/genome_datasets/Olea_europaea_all_scaffolds.fa.gz 
# 	gunzip -c {output.gz} > {output.uncmp}
# 	'''
# 
# rule download_busco_db:
#     output: 
#         db='inputs/busco/embryophyta_odb9/',
#         tgz='inputs/busco/embryophyta_odb9.tar.gz'
#     shell:'''
#     wget -O {output.tgz} http://busco.ezlab.org/datasets/embryophyta_odb9.tar.gz
# 	mkdir -p {output.db}
#     tar -xf {output.tgz} --strip-components=1 -C {output.db}
#     '''
# rule run_busco:
#     output: 'outputs/busco/run_wild_olive_busco'
#     input: 
#         genome='inputs/sylvestris/Olea_europaea_all_scaffolds.fa',
#         busco_db='inputs/busco/embryophyta_odb9'
#     conda:  "envs/busco.yml"
#     shell:'''
# 	run_busco -i {input.genome} -o wild_olive_busco -l {input.busco_db} -m geno
#         mv run_wild_olive_busco {output}
# '''