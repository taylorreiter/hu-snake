import pysam

ENV = "clustering-env.yml"
PFAM_LINK = "https://pfam.xfam.org/family/PF00521/alignment/full" # link to the PFAM seqs to use; here gyrA
PFAM_BASE = "PF00521_full_gyra" # file basename to save PFAM as
FAA = "loosetrim_plus_bin" # input basename faa amino acid sequence
OUT_BASE = "loosetrim-plus-bin-PF00521-hmmscanT100" # input basename for hmmer output

rule all:
    input: 
        expand("outputs/pid/{out_base}.mat", out_base = OUT_BASE)
 
rule download_pfam:
    output: "inputs/pfam/{pfam_base}.sto"
    params: 
        pfam_link = PFAM_LINK,
        pfam_base = PFAM_BASE,
        out_dir = "inputs/pfam"
    shell:'''
    wget -O {params.out_dir}/{params.pfam_base}.sto {params.pfam_link}
    '''

rule hmmbuild:
    output: "outputs/hmmbuild/{pfam_base}.hmm"
    input: "inputs/pfam/{pfam_base}.sto"
    conda: ENV
    shell:'''
    hmmbuild {output} {input}
    hmmpress {output}
    '''

rule download_faa:
    output: "inputs/plass/{faa}.faa"
    params: faa = FAA
    shell:'''
    # add download link from osf
    #mv {params.faa}.faa {output}
    touch {output}
    '''

rule format_faa_headers_cut:
    output: "outputs/plass/{faa}.cut.faa"
    input: "inputs/plass/{faa}.faa"
    shell:'''
    cut -d ' ' -f1 {input} > {output}
    '''

rule format_faa_headers_dup:
    output: "outputs/plass/{faa}.cut.dup.faa"
    input: "outputs/plass/{faa}.cut.faa"
    shell:'''
    awk '(/^>/ && s[$0]++){{$0=$0"_"s[$0]}}1;' {input} > {output}
    '''
 
rule hmmscan:
    output:
        out = "outputs/hmmscan/{out_base}.out",
        tbl = "outputs/hmmscan/{out_base}-tbl.out",
        dom = "outputs/hmmscan/{out_base}-dom.out"
    input:
        hmm = expand("outputs/hmmbuild/{pfam_base}.hmm", pfam_base = PFAM_BASE),
        faa = expand("outputs/plass/{faa}.cut.dup.faa", faa = FAA)
    conda: ENV
    shell:'''
    hmmscan -T 100 -o {output.out} --tblout {output.tbl} --domtblout {output.dom} {input.hmm} {input.faa} 
    '''

rule get_rhmmer:
    output: "rhmmer.R"
    shell:'''
    wget -O rhmmer.R https://raw.githubusercontent.com/arendsee/rhmmer/master/R/parse.R
    '''

rule def_overlap_window_and_find_reads:
    input:
        dom = "outputs/hmmscan/{out_base}-dom.out",
        rhmmer = "rhmmer.R"
    output:
        keep = "outputs/hmmscan/{out_base}-NAMES.txt"
    conda: ENV
    script:'clustering-workflow-overlaps.R'

rule grab_overlap_faa:
    output: filtered_fa = "outputs/hmmscan/{out_base}.faa"
    input:
        fasta = expand("outputs/plass/{faa}.cut.dup.faa", faa = FAA),
        names = "outputs/hmmscan/{out_base}-NAMES.txt" 
    shell:'''
    ./extract-hmmscan-matches.py {input.names} {input.fasta} > {output}
    '''

rule mafft:
    output: "outputs/mafft/{out_base}-mafft.faa"
    input: "outputs/hmmscan/{out_base}.faa"
    conda: ENV
    shell:'''
    mafft --auto --reorder {input} > {output}
    '''

rule mafft_to_sto:
    output: "outputs/mafft/{out_base}.sto"
    input: "outputs/mafft/{out_base}-mafft.faa"
    run:
        from Bio import SeqIO
        from Bio import AlignIO

        align = AlignIO.read(str(input), "fasta")
        with open(str(output), "w") as handle:
            count = SeqIO.write(align, handle, "stockholm")

rule calc_pid:
    output: "outputs/pid/{out_base}.pid"
    input: "outputs/mafft/{out_base}.sto"
    conda: ENV
    shell:'''
    esl-alipid {input} > {output}
    '''

rule pid_mat:
    output:
        mat = "outputs/pid/{out_base}.mat" 
    input:
        pid = "outputs/pid/{out_base}.pid"
    conda: ENV
    script: 'clustering-workflow-mat.R'

#rule pid_plot:
#    output:
#        html = 
#        pdf =
#    input: 
#        mat = "outputs/pid/{out_base}.mat"
#    conda: ENV
#    script: 'clustering-workflow-plot.R'



#    run:
#        # inspired by https://gist.github.com/ngcrawford/3217022
#        def splitIterator(text, size):
#            """
#            Iterator that splits string into list of substrings.
#            """
#            assert size > 0, "size should be > 0"
#            for start in range(0, len(text), size):
#                yield text[start:start + size]
#        
#
#        # Filter names for contings
#        names = open(input.names, "r")
#        seq_names = names.read().split('\n')
#        names.close()
#
#        # Write contigs/chrms to output fasta
#        fasta = pysam.Fastafile(str(input.fasta))
#        filtered_fa = open(output.filtered_fa,'w')
#  
#        for name in seq_names:
#            seq = fasta.fetch(name)
#            filtered_fa.write(">" + name + "\n")
#            # split lines on 80 characters
#            [filtered_fa.write(chars + "\n" ) for chars in splitIterator(seq, 80)]
#
#        # Clean up open files
#        filtered_fa.close() 
