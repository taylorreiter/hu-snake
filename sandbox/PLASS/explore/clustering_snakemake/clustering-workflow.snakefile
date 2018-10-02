import pysam

# upstream code:
# cp ../../../../outputs/hu-bins/prokka/sb1_tmp/all_sb1_bin_prokka.faa .
# sed 's/>/&BIN-/' all_sb1_bin_prokka.faa > all_sb1_bin_prokka2.faa
# cat all_hardtrim.plass.c100.faa all_sb1_bin_prokka2.faa > all_hardtrim.plass.c100.all_bin.faa
# cat all_loosetrim.plass.c100.faa all_sb1_bin_prokka2.faa > all_loosetrim.plass.c100.all_bin.faa

ENV = "clustering-env.yml"

PFAM_LINK = "https://pfam.xfam.org/family/PF00521/alignment/full" # link to the PFAM seqs to use; here gyrA
PFAM_BASE = "PF00521_full_gyra" # file basename to save PFAM as

## original gyrA:
#FAA = "loosetrim_plus_bin" # input basename faa amino acid sequence
#OUT_BASE = "loosetrim-plus-bin-PF00521-hmmscanT100" # input basename for hmmer output

## new gyrA (has full loose trim & hard trim, as well as concatenated seqs from bins)
FAA = "all_hardtrim.plass.c100.all_bin"
OUT_BASE = "plass-hardtrim-all-bin-PF00521-hmmscanT100"
#FAA = "all_loosetrim.plass.c100.all_bin"
#OUT_BASE = "plass-loosetrim-all-bin-PF00521-hmmscanT100"

#PFAM_LINK = "https://pfam.xfam.org/family/PF00204/alignment/full" #gyrB
#PFAM_BASE = "PF00204_full_gyrb"
#FAA = "all_hardtrim.plass.c100"
#OUT_BASE = "plass-hardtrim-PF00204-hmmscanT100"
#FAA = "all_loosetrim.plass.c100"
#OUT_BASE = "plass-loosetrim-PF00204-hmmscanT100"

#PFAM_LINK = "https://pfam.xfam.org/family/PF00181/alignment/full" # rplb
#PFAM_BASE = "PF00181_full_rplb"
#FAA = "all_hardtrim.plass.c100"
#OUT_BASE = "plass-hardtrim-PF00181-hmmscanT100"
#FAA = "all_loosetrim.plass.c100"
#OUT_BASE = "plass-loosetrim-PF00181-hmmscanT100"

#PFAM_LINK = "https://pfam.xfam.org/family/PF00189/alignment/full" # rpsc
#PFAM_BASE = "PF00189_full_rpsc"
#FAA = "all_hardtrim.plass.c100"
#OUT_BASE = "plass-hardtrim-PF00189-hmmscanT100"
#FAA = "all_loosetrim.plass.c100"
#OUT_BASE = "plass-loosetrim-PF00189-hmmscanT100"

#PFAM_LINK = "https://pfam.xfam.org/family/PF00154/alignment/full" # recA
#PFAM_BASE = "PF00154_full_reca"
#FAA = "all_hardtrim.plass.c100"
#OUT_BASE = "plass-hardtrim-PF00154-hmmscanT100"
#FAA = "all_loosetrim.plass.c100"
#OUT_BASE = "plass-loosetrim-PF00154-hmmscanT100"

#PFAM_LINK = "https://pfam.xfam.org/family/PF01411/alignment/full" # alaS
#PFAM_BASE = "PF01411_full_alas"
#FAA = "all_hardtrim.plass.c100"
#OUT_BASE = "plass-hardtrim-PF01411-hmmscanT100"
#FAA = "all_loosetrim.plass.c100"
#OUT_BASE = "plass-loosetrim-PF01411-hmmscanT100"

#PFAM_LINK = "https://pfam.xfam.org/family/PF00562/alignment/full" # rpb2 d6
#PFAM_BASE = "PF00562_full_rpb2d6"
#FAA = "all_hardtrim.plass.c100"
#OUT_BASE = "plass-hardtrim-PF00562-hmmscanT100"
#FAA = "all_loosetrim.plass.c100"
#OUT_BASE = "plass-loosetrim-PF00562-hmmscanT100"

rule all:
    input: 
        expand("outputs/pid/{out_base}-mds.csv", out_base = OUT_BASE)
 
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

rule pid_mds:
    output:
        mds = "outputs/pid/{out_base}-mds.csv" 
    input:
        mat = "outputs/pid/{out_base}.mat", 
        info = "inputs/hu_info_sb1.csv"
    conda: ENV
    script: 'clustering-workflow-mds.R'
