# upstream code:
# cp ../../../../outputs/hu-bins/prokka/sb1_tmp/all_sb1_bin_prokka.faa .
# sed 's/>/&BIN-/' all_sb1_bin_prokka.faa > all_sb1_bin_prokka2.faa
# cat all_hardtrim.plass.c100.faa all_sb1_bin_prokka2.faa > all_hardtrim.plass.c100.all_bin.faa
# cat all_loosetrim.plass.c100.faa all_sb1_bin_prokka2.faa > all_loosetrim.plass.c100.all_bin.faa

ENV = "clustering-env.yml"

PFAM_BASE=["PF00521_gyra",
            "PF00204_gyrb",
            "PF00181_rplb", 
            "PF00189_rpsc",
            "PF00154_reca",
            "PF01411_alas",
            "PF00562_rpb2d6"] 

FAA = ["all_hardtrim.plass.c100.all_bin"]
#FAA = "all_loosetrim.plass.c100.all_bin"

TRIM = ["hard"]
# TRIM = ["loose"]

rule all:
    input: 
        expand("outputs/pid/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-mds.csv", trim = TRIM, pfam_base = PFAM_BASE)

rule download_pfam:
    output: "inputs/pfam/{pfam_base}.sto"
    params: 
        out_dir = "inputs/pfam"
    shell:'''
    pfam=$(echo {wildcards.pfam_base} | cut -f1 -d"_")
    wget -O {output} https://pfam.xfam.org/family/${{pfam}}/alignment/full
    '''

rule hmmbuild:
    input: "inputs/pfam/{pfam_base}.sto"
    output: "outputs/hmmbuild/{pfam_base}.hmm"
    shell: '''
    hmmbuild outputs/hmmbuild/{wildcards.pfam_base}.hmm inputs/pfam/{wildcards.pfam_base}.sto   
    hmmpress outputs/hmmbuild/{wildcards.pfam_base}.hmm
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
        out = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.out",
        tbl = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-tbl.out",
        dom = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-dom.out"
    input:
        hmm = "outputs/hmmbuild/{pfam_base}.hmm",
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
        dom = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-dom.out",
        rhmmer = "rhmmer.R"
    output:
        keep = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-NAMES.txt"
    conda: ENV
    script:'clustering-workflow-overlaps.R'

rule grab_overlap_faa:
    output: filtered_fa = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.faa"
    input:
        fasta = expand("outputs/plass/{faa}.cut.dup.faa", faa = FAA),
        names = "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-NAMES.txt" 
    shell:'''
    ./extract-hmmscan-matches.py {input.names} {input.fasta} > {output}
    '''

rule mafft:
    output: "outputs/mafft/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-mafft.faa"
    input: "outputs/hmmscan/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.faa"
    conda: ENV
    shell:'''
    mafft --auto --reorder {input} > {output}
    '''

rule mafft_to_sto:
    output: "outputs/mafft/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.sto"
    input: "outputs/mafft/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-mafft.faa"
    run:
        from Bio import SeqIO
        from Bio import AlignIO

        align = AlignIO.read(str(input), "fasta")
        with open(str(output), "w") as handle:
            count = SeqIO.write(align, handle, "stockholm")

rule calc_pid:
    output: "outputs/pid/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.pid"
    input: "outputs/mafft/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.sto"
    conda: ENV
    shell:'''
    esl-alipid {input} > {output}
    '''

rule pid_mat:
    output:
        mat = "outputs/pid/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.mat" 
    input:
        pid = "outputs/pid/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.pid"
    conda: ENV
    script: 'clustering-workflow-mat.R'

rule pid_mds:
    output:
        mds = "outputs/pid/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100-mds.csv" 
    input:
        mat = "outputs/pid/plass-{trim}trim-all-bin-{pfam_base}-hmmscanT100.mat", 
        info = "inputs/hu_info_sb1.csv"
    conda: ENV
    script: 'clustering-workflow-mds.R'
