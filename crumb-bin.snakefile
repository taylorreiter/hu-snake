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
    curl -L -o {output} https://osf.io/u5yqf/download
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

rule download_blast_db_seqs_tax:
    output: 'inputs/blast_db/taxonomy-aa.faa'
    shell:'''
    curl -L -o gyrA.faa 'http://www.uniprot.org/uniprot/?query=gyra&sort=score&format=fasta'
    curl -L -o gyrB.faa 'http://www.uniprot.org/uniprot/?query=gyrb&sort=score&format=fasta'
    curl -L -o recA.faa 'http://www.uniprot.org/uniprot/?query=reca&sort=score&format=fasta'
    cat gyrA.faa gyrB.faa recA.faa > {output}
    rm gyrA.faa gyrB.faa recA.faa
    '''
    
rule make_blast_db_tax:
    output: 'inputs/blast_db/taxonomy-aa.faa.psq'
    input: 'inputs/blast_db/taxonomy-aa.faa'
    conda: 'env.yml'
    shell:'''
    makeblastdb -in {input} -dbtype prot
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
    Rscript --vanilla merge_fasta.R {input.mh} {input.uni} {output}
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

# busco ------------------------------------------------------------

# estimate number of single copy orthologs that are contained in the unitigs, subtracts, and in the assemblies originating from the crumbs_bin.
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
    output: 'outputs/hu-crumbs_bin/unitigs/busco/run_{crumb_bin}_bac'
    input: 
        crumb_bin_in='inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.contigs.fa.gz.crumb_bin.fa',
        busco_db='inputs/busco/bacteria_odb9/'
    conda:  "env.yml"
    shell:'''
	run_busco -i {input.crumb_bin_in} -o {wildcards.crumb_bin}_bac -l {input.busco_db} -m geno
    mv run_{wildcards.crumb_bin}_bac {output}
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
    Rscript --vanilla merge_fasta.R {input.mh} {input.uni} {output}
    '''

# blast ------------------------------------------------------------------

rule blastp_subtracts:
    output: 'outputs/hu-crumbs_bin/subtracts/blast/{crumb_bin}-blastp.tab'
    input: 
        db = 'inputs/blast_db/interesting-aa.faa.psq',
        query = 'outputs/hu-crumbs_bin/subtracts/prokka-all/{crumb_bin}.faa'
    conda: 'env.yml'   
    shell: '''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
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
    
# blast ------------------------------------------------------------------

# rule blastp_assemblies:
#     output: 'outputs/hu-crumbs_bin/assembly/blast/{crumb_bin}-blastp.tab'
#     input: 
#         db = 'inputs/blast_db/interesting-aa.faa.psq',
#         query = 'outputs/hu-crumbs_bin/assembly/prokka/{crumb_bin}.faa'
#     conda: 'env.yml'   
#     shell: '''
#     touch {input.db}
#     blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -out {output} -outfmt 6
#     '''

rule blastp_tax_assemblies:
    output: 'outputs/hu-crumbs_bin/assembly/blast/{crumb_bin}-tax-blastp.tab'
    input:
        db='inputs/blast_db/taxonomy-aa.faa.psq',
        query='outputs/hu-crumbs_bin/assembly/prokka/{crumb_bin}.faa'
    conda: 'env.yml'
    shell:'''
    touch {input.db}
    blastp -query {input.query} -db inputs/blast_db/interesting-aa.faa -evalue 1E-10 -max_target_seqs 5 -out {output} -outfmt 6
    '''
# busco ------------------------------------------------------------------

rule run_busco_bac_assembly:
    output: 'outputs/hu-crumbs_bin/assembly/busco/run_{crumb_bin}_bac'
    input: 
        crumb_bin_in='inputs/hu-crumbs_bin/{crumb_bin}.fa.cdbg_ids.reads.fa.gz.crumb_bin.fa.assembly.fa',
        busco_db='inputs/busco/bacteria_odb9/'
    conda:  "env.yml"
    shell:'''
	run_busco -i {input.crumb_bin_in} -o {wildcards.crumb_bin}_bac -l {input.busco_db} -m geno
    mv run_{wildcards.crumb_bin}_bac {output}
    '''
