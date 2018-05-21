from snakemake.utils import R
configfile: "config.yml"

# 1. 0.fa.cdbg_ids.contigs.fa.gz.croissant.fa
# 2. 0.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa
# 3. 0.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa

# 1 is the unitigs from the croissant
# 2 is #1 - #3
# 3 is the assembly of the croissant reads using megahit

rule download_croissants:
    output: 'inputs/hu-croissants/hu-croissants.tar.gz'
    shell:'''
    curl -L -o {output} https://osf.io/u5yqf/download
    '''

rule unpack_croissants:
    output: 
        dynamic('inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa'),
        dynamic('inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa'),
        dynamic('inputs/hu-croissants/{croissant}.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa')
    input: 'inputs/hu-croissants/hu-croissants.tar.gz'
    params: output_folder = 'inputs/hu-croissants/'
    shell:'''
    tar xf {input} --directory {params.output_folder}
    '''

# resources for all analyses:

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
    
# SUMMARY STATS ----------------------------------------------------------

rule index_hu_unitigs:
    output: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.fai'
    input: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
    
rule index_hu_assembly:
    output: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa.fai'
    input: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
   
rule index_hu_subtract:
    output: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa.fai'
    input: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa'
    conda: 'env.yml'
    shell:'''
    samtools faidx {input}
    '''
rule summarize_hu:
    output: 'outputs/hu-croissants/summary_of_inputs.tsv'
    input: 
        dynamic('inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.fai'),
        dynamic('inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa.fai'),
        dynamic('inputs/hu-croissants/{croissant}.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa.fai')
    conda: 'env2.yml'
    shell:'''
    Rscript --vanilla scripts/skim_input.R
    '''

# UNITIGS ################################################################   
# megahit & annotate unitigs ---------------------------------------------
    
rule assemble_croissant_unitigs:
    output: 'outputs/hu-croissants/unitigs/megahit/{croissant}.contigs.fa'
    input: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/unitigs/megahit'
    shell:'''
    # megahit does not allow force overwrite, so each assembly needs to take place in it's own directory.
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.croissant} --out-prefix {wildcards.croissant} 
    # move the final assembly to a folder containing all assemblies
    mv {wildcards.croissant}/{wildcards.croissant}.contigs.fa {params.output_folder}/{wildcards.croissant}.contigs.fa
    # remove the original megahit assembly folder, which is in the main directory.
    rm -rf {wildcards.croissant}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_croissant_unitigs:
    output: 'outputs/hu-croissants/unitigs/megahit-prokka/{croissant}.faa'
    input:  'outputs/hu-croissants/unitigs/megahit/{croissant}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/unitigs/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.croissant} --metagenome --force
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------

rule prokka_croissants_unitigs:
    output: 'outputs/hu-croissants/unitigs/unitig-prokka/{croissant}.faa'
    input:  'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/unitigs/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.croissant} --metagenome --force
    touch {output}
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_croissant_unitigs:
    output: 'outputs/hu-croissants/unitigs/prokka-all/{croissant}.faa'
    input:  
        mh = 'outputs/hu-croissants/unitigs/megahit-prokka/{croissant}.faa',
        uni = 'outputs/hu-croissants/unitigs/unitig-prokka/{croissant}.faa'
    shell:'''
    Rscript --vanilla merge_fasta.R {input.mh} {input.uni} {output}
    '''

# blast ------------------------------------------------------------------

rule blastp_croissant_unitigs:
    output: 'outputs/hu-croissants/unitigs/blast/{croissant}-blastp.tab'
    input: 
        db = 'inputs/blast_db/interesting-aa.faa.psq',
        query = 'outputs/hu-croissants/unitigs/prokka-all/{croissant}.faa'
    conda: 'env.yml'   
    shell: '''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
    '''

# busco ------------------------------------------------------------

# estimate number of single copy orthologs that are contained in the unitigs, subtracts, and in the assemblies originating from the croissants.
# # https://www.microbe.net/2017/12/13/why-genome-completeness-and-contamination-estimates-are-more-complicated-than-you-think/

rule download_busco_bac:
    output: 
        db='inputs/busco/bacteria_odb9/',
        tgz='inputs/busco/bacteria_odb9.tar.gz'
    shell:'''
    wget -O {output.tgz} http://busco.ezlab.org/datasets/bacteria_odb9.tar.gz
	mkdir -p {output.db}
    tar -xf {output.tgz} --strip-components=1 -C {output.db}
    '''

rule run_busco_bac_unitigs:
    output: 'outputs/hu-croissants/unitigs/busco/run_{croissant}_bac'
    input: 
        croissant_in='inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa',
        busco_db='inputs/busco/bacteria_odb9/'
    conda:  "env.yml"
    shell:'''
	run_busco -i {input.croissant_in} -o {wildcards.croissant}_bac -l {input.busco_db} -m geno
    mv run_{wildcards.croissant}_bac {output}
    '''
    
# SUBTRACTION ################################################################

rule assemble_croissant_subtracts:
    output: 'outputs/hu-croissants/subtracts/megahit/{croissant}.contigs.fa'
    input: 'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/subtracts/megahit'
    shell:'''
    megahit -r {input} --min-contig-len 142 --out-dir {wildcards.croissant} --out-prefix {wildcards.croissant} 
    mv {wildcards.croissant}/{wildcards.croissant}.contigs.fa {params.output_folder}/{wildcards.croissant}.contigs.fa
    rm -rf {wildcards.croissant}
    ''' 

# annotate megahit assemblies and hu genomes with prokka
rule prokka_megahit_croissants_subtracts:
    output: 'outputs/hu-croissants/subtracts/megahit-prokka/{croissant}.faa'
    input:  'outputs/hu-croissants/subtracts/megahit/{croissant}.contigs.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/subtracts/megahit-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.croissant} --metagenome --force
    touch {output}
    '''

# annotate unitigs ------------------------------------------------------
    
rule prokka_unitig_croissants_subtracts:
    output: 'outputs/hu-croissants/subtracts/unitig-prokka/{croissant}.faa'
    input:  'inputs/hu-croissants/{croissant}.fa.cdbg_ids.contigs.fa.gz.croissant.fa.sub.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/subtracts/unitig-prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.croissant} --metagenome --force
    touch {output}
    '''
       
# combine megahit and unitig annotations ---------------------------------

rule combine_prokka_subtracts:
    output: 'outputs/hu-croissants/subtracts/prokka-all/{croissant}.faa'
    input:  
        mh = 'outputs/hu-croissants/subtracts/megahit-prokka/{croissant}.faa',
        uni = 'outputs/hu-croissants/subtracts/unitig-prokka/{croissant}.faa'
    shell:'''
    Rscript --vanilla merge_fasta.R {input.mh} {input.uni} {output}
    '''

# blast ------------------------------------------------------------------

rule blastp_subtracts:
    output: 'outputs/hu-croissants/subtracts/blast/{croissant}-blastp.tab'
    input: 
        db = 'inputs/blast_db/interesting-aa.faa.psq',
        query = 'outputs/hu-croissants/subtracts/prokka-all/{croissant}.faa'
    conda: 'env.yml'   
    shell: '''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
    '''

# ASSEMBLIES ################################################################

# annotate megahit assemblies and hu genomes with prokka
rule prokka_croissants_assemblies:
    output: 'outputs/hu-croissants/assembly/prokka/{croissant}.faa'
    input:  'inputs/hu-croissants/{croissant}.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa'
    conda: 'env.yml'
    params:
        output_folder = 'outputs/hu-croissants/assembly/prokka'
    shell:'''
    prokka {input} --outdir {params.output_folder} --prefix {wildcards.croissant} --metagenome --force
    touch {output}
    '''
    
# blast ------------------------------------------------------------------

rule blastp_assemblies:
    output: 'outputs/hu-croissants/assembly/blast/{croissant}-blastp.tab'
    input: 
        db = 'inputs/blast_db/interesting-aa.faa.psq',
        query = 'outputs/hu-croissants/assembly/prokka/{croissant}.faa'
    conda: 'env.yml'   
    shell: '''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
    '''
# busco ------------------------------------------------------------------

rule run_busco_bac_assembly:
    output: 'outputs/hu-croissants/assembly/busco/run_{croissant}_bac'
    input: 
        croissant_in='inputs/hu-croissants/{croissant}.fa.cdbg_ids.reads.fa.gz.croissant.fa.assembly.fa',
        busco_db='inputs/busco/bacteria_odb9/'
    conda:  "env.yml"
    shell:'''
	run_busco -i {input.croissant_in} -o {wildcards.croissant}_bac -l {input.busco_db} -m geno
    mv run_{wildcards.croissant}_bac {output}
    '''
