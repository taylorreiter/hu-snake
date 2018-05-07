from snakemake.utils import R
configfile: "config.yml"

rule download_croissants:
    output: 'inputs/hu-croissants/{sample}.fa.nbhd.fa.gz'
    shell:'''
    curl -L -o inputs/hu-croissants/hu-croissant08.fa.nbhd.fa.gz https://osf.io/g2zvh/download 
    curl -L -o inputs/hu-croissants/hu-croissant11.fa.nbhd.fa.gz https://osf.io/5k2zj/download
    curl -L -o inputs/hu-croissants/hu-croissant19.fa.nbhd.fa.gz https://osf.io/29bk8/download
    ''' 
    
# MEGAHIT & ANNOTATE UNITIGS ---------------------------------------------
    
rule assemble_croissants:
    output: 'outputs/hu-croissants/megahit/{sample}.contigs.fa'
    input: 'inputs/hu-croissants/{sample}.fa.nbhd.fa.gz'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 141 --out-dir {wildcards.sample} --out-prefix {wildcards.sample} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.sample}/{wildcards.sample}.contigs.fa {params.output_folder}/{wildcards.sample}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.sample}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_croissants:
    output: 'outputs/hu-croissants/megahit-prokka/{sample}.faa'
    input:  'outputs/hu-croissants/megahit/{sample}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.sample} --metagenome --force
    touch {output}
    '''

# ANNOTATE UNITIGS ------------------------------------------------------

rule gunzip_unitig_croissants:
    output: 'inputs/hu-croissants/{sample}.fa.nbhd.fa'
    input: 'inputs/hu-croissants/{sample}.fa.nbhd.fa.gz'
    shell:'''
    gunzip {input}
    '''
    
rule prokka_unitig_croissants:
    output: 'outputs/hu-croissants/unitig-prokka/{sample}.faa'
    input:  'inputs/hu-croissants/{sample}.fa.nbhd.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.sample} --metagenome --force
    touch {output}
    '''
       
# COMBINE MEGAHIT & UNITIG ANNOTATIONS ----------------------------------

rule combine_prokka:
    output: 'outputs/hu-croissants/prokka-all/{sample}.faa'
    input:  
        mh = 'outputs/hu-croissants/megahit-prokka/{sample}.faa',
        uni = 'outputs/hu-croissants/unitig-prokka/{sample}.faa'
    shell:'''
    Rscript --vanilla merge_fasta.R {input.mh} {input.uni} {output}
    '''

# BLAST -----------------------------------------------------------------

# blast the amino acid sequences against interesting dbs

rule download_blast_db_seqs:
    output: 'inputs/blast_db/interesting-aa.faa'
    shell:'''
    # Hydrocarbon degradation
    curl -L -o alkylsuccinate-synthase.faa 'http://www.uniprot.org/uniprot/?query=alkylsuccinate+synthase&sort=score&format=fasta'
    curl -L -o benzylsuccinate-synthase.faa 'http://www.uniprot.org/uniprot/?query=benzylsuccinate+synthase&sort=score&format=fasta'
    curl -L -o succinyl-CoA-R-Benzylsuccinate-CoA-transferase.faa 'http://www.uniprot.org/uniprot/?query=succinyl-CoA+Benzylsuccinate-CoA+transferase&sort=score&format=fasta'
    curl -L -o benzylsuccinyl-CoA-dehydrogenase.faa 'http://www.uniprot.org/uniprot/?query=benzylsuccinyl-CoA+dehydrogenase&sort=score&format=fasta'
    curl -L -o benzoyl-CoA-reductase.faa 'http://www.uniprot.org/uniprot/?query=benzoyl-CoA-reductase&sort=score&format=fasta'
    # original paper had an 'e' on benzoyl; I removed this. 
    curl -L -o alkane-1-monooxygenase.faa 'http://www.uniprot.org/uniprot/?query=alkane-1-monooxygenase&sort=score&format=fasta'
    curl -L -o benzene-carboxylase.faa 'http://www.uniprot.org/uniprot/?query=benzene+carboxylase&sort=score&format=fasta'
    curl -L -o monooxygenase-flavin-binding-protein.faa 'http://www.uniprot.org/uniprot/?query=monooxygenase+flavin-binding+protein&sort=score&format=fasta'
    curl -L -o ethylbenzene-dehydrogenase.faa 'http://www.uniprot.org/uniprot/?query=ethylbenzene+dehydrogenase&sort=score&format=fasta'
    curl -L -o cyclohex-1-5-diene-1-carboxyl-CoA-hydratase.faa 'http://www.uniprot.org/uniprot/?query=cyclohex-1,5-diene-1-carboxyl-CoA+hydratase&sort=score&format=fasta'
    curl -L -o 6-hydroxycyclohex-1-ene-1-carboxyl-CoA-dehydrogenase.faa 'http://www.uniprot.org/uniprot/?query=6-hydroxycyclohex-1-ene-1-carboxyl-CoA+dehydrogenase&sort=score&format=fasta'
    curl -L -o benzoate-CoA-ligase.faa 'http://www.uniprot.org/uniprot/?query=benzoate-CoA+ligase&sort=score&format=fasta'
    curl -L -o benzoate-CoA-reductase.faa 'http://www.uniprot.org/uniprot/?query=benzoate-CoA+reductase&sort=score&format=fasta'
    curl -L -o cyclohex-1-ene-1-carboxyl-CoA-hydrolase.faa 'http://www.uniprot.org/uniprot/?query=cyclohex-1-ene-1-carboxyl-CoA+hydrolase&sort=score&format=fasta'
    curl -L -o 2-hydroxycyclohexane-carboxyl-CoA-dehydrogenase.faa 'http://www.uniprot.org/uniprot/?query=2-hydroxycyclohexane+carboxyl-CoA+dehydrogenase&sort=score&format=fasta'
    curl -L -o 2-ketocyclohexane-carboxyl-CoA-hydrolase.faa 'http://www.uniprot.org/uniprot/?query=2-ketocyclohexane+carboxyl-CoA+hydrolase&sort=score&format=fasta'
    # Sulfite reduction
    curl -L -o sulfate-adenylyl-transferase.faa 'http://www.uniprot.org/uniprot/?query=sulfate+adenylyl+transferase&sort=score&format=fasta'
    curl -L -o sulfite-reductase-dissimilatory.faa 'http://www.uniprot.org/uniprot/?query=sulfite+reductase+dissimilatory&sort=score&format=fasta'
    # AromaDeg
    curl -L -o AromaDeg_all_unaligned.fasta.bz2 http://aromadeg.siona.helmholtz-hzi.de/database/AromaDeg_all_unaligned.fasta.bz2
    bzip2 -d AromaDeg_all_unaligned.fasta.bz2
    cat *faa *fasta > {output}
    rm *faa *fasta*
    '''

rule make_blast_db:
    output: 'inputs/blast_db/interesting-aa.faa.psq'
    input: 'inputs/blast_db/interesting-aa.faa'
    conda: 'env.yml'
    shell:'''
    makeblastdb -in {input} -dbtype prot
    '''

rule blastp:
    output: 'outputs/hu-croissants/blast/{sample}-blastp.tab'
    input: 
        db = 'inputs/blast_db/interesting-aa.faa.psq',
        query = 'outputs/hu-croissants/prokka-all/{sample}.faa'
    conda: 'env.yml'   
    shell: '''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
    '''

# BUSCO ------------------------------------------------------------

# NOTE this relies on the dbs downloaded in genome.snakefile. Fix this later so input files will be found.

rule run_busco_bac_c:
    output: 'outputs/hu-croissants/busco/run_{sample}_bac'
    input: 
        croissant_in='inputs/hu-croissants/{sample}.fa.nbhd.fa',
        busco_db='inputs/busco/bacteria_odb9/'
    conda:  "env.yml"
    shell:'''
	run_busco -i {input.croissant_in} -o {wildcards.sample}_bac -l {input.busco_db} -m geno
    mv run_{wildcards.sample}_bac {output}
    '''

# rule run_busco_arch:
#     output: 'outputs/hu-croissants/busco/run_{croissant}_arch'
#     input: 
#         croissant_in='inputs/hu-croissants/{croissant}.fa',
#         busco_db='inputs/busco/archea_odb9/'
#     conda:  "env.yml"
#     shell:'''
# 	run_busco -i {input.croissant_in} -o {wildcards.croissant}_arch -l {input.busco_db} -m geno
#     mv run_{wildcards.croissant}_arch {output}
#     ''' 