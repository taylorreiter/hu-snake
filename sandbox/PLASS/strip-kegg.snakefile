ENV = "env-kegg.yml"

rule all:
    input:
        'outputs/GhostKOALA/other/crumb_nitro_genes.csv',
        dynamic('outputs/GhostKOALA-plots/{crumb_bin}.pdf'),
        #'outputs/GhostKOALA-plots/kegg-each.pdf',
        #'outputs/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        #'outputs/GhostKOALA-plots/kegg-all-col-kingdom.pdf',
        'outputs/GhostKOALA/an-et-al/an_ass.tsv',
        'outputs/GhostKOALA/an-et-al/an_ass.pdf',
        #'outputs/GhostKOALA/marker-genes/blast_aa.csv',
        #'outputs/GhostKOALA/marker-genes/blast_aas_full.csv',
        'outputs/GhostKOALA/marker-genes/novel_crumb_marks.csv',
        'outputs/GhostKOALA/other/crumb_nitro_genes.csv',
        'outputs/figures/fig5a.pdf',
        'outputs/figures/fig5a.png',
        'outputs/figures/fig5b.pdf',
        'outputs/figures/fig5b.png'

# CRUMB ASSEMBLIES #################################################

rule download_kegg:
    output: 
        ann = 'outputs/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/GhostKOALA/user.out.top'  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

rule split_kegg:
    output:
        ann = dynamic('outputs/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        tax = dynamic('outputs/GhostKOALA/{crumb_bin}.ko-tax')  
    input: 
        ann = 'outputs/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/GhostKOALA/user.out.top'  
    run:
        import pandas as pd
        hu = pd.read_table(str(input.ann), header=None)
        hu[0] = hu[0].replace("user:",  "", regex = True)
        hu[6] = hu[0].replace("_SRR1976948.[0-9_]{2,12}", "", regex = True)
        for sample, hu_bin in hu.groupby(6):
            hu_bin.to_csv(f"outputs/GhostKOALA/{sample}.ko-ann-full.txt", index = False, header = False, sep = "\t")

        tax = pd.read_table(str(input.tax), header=None)
        tax[0] = tax[0].replace("user:",  "", regex = True)
        tax[7] = tax[0].replace("_SRR1976948.[0-9_]{2,12}", "", regex = True)
        for sample, tax_bin in tax.groupby(7):
            tax_bin.to_csv(f"outputs/GhostKOALA/{sample}.ko-tax", index = False, header = False, sep = "\t")

# PLOTS

rule plot_kegg_taxonomy:
    output: 'outputs/GhostKOALA-plots/{crumb_bin}.pdf'
    input: 'outputs/GhostKOALA/{crumb_bin}.ko-tax'
    run:
        import pandas as pd
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.rcParams['figure.figsize'] = [17, 8]
        
        hu = pd.read_table(str(input), header=None)
        hu = hu.loc[hu[6] > 100] # recommended by GhostKOALA
        hu = hu[~hu[4].str.contains("candidate")]
        hu = hu[~hu[4].str.contains("Candidatus")]
        
        tax = hu.groupby(4)
        tax = tax.filter(lambda x: len(x) > (hu.shape[0]/1500)) # dynamically adjust values on x axis
        
        plot = tax[4].value_counts().plot(kind = 'bar')
        plt.tight_layout()
        plot.get_figure().savefig(str(output), format='pdf') 

rule plot_kegg_annotations:
    output: 
        plot1='outputs/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        plot2b= 'outputs/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        info = 'inputs/hu_info.csv'
    params: dir = "outputs/GhostKOALA"
    conda: ENV
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations:
    output:
        tsv='outputs/GhostKOALA/an-et-al/an_ass.tsv',
        plot='outputs/GhostKOALA/an-et-al/an_ass.pdf'
    input:
        an ='../an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/GhostKOALA/{crumb_bin}.ko-ann-full.txt')
    params: 
        title = "Assembly",
        dir = 'outputs/GhostKOALA/'
    conda: ENV
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''

# HU BINS ############################################

rule download_kegg_bins:
    output: 
        ann = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-bins/GhostKOALA/user.out.top'  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

# COMBINED ##########################################

rule blast_marker_genes:
    output: 
        blast_aas = 'outputs/GhostKOALA/marker-genes/blast_aa.csv',
        blast_aas_full = 'outputs/GhostKOALA/marker-genes/blast_aas_full.csv'
    input:
        info = 'inputs/hu_info.csv',
        crumb_kegg = 'outputs/GhostKOALA/user_ko_definition.txt',
        bin_kegg = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt',
        faa = 'outputs/hu-genomes-plass-clean/sb1-plass.nbhd.nostop.99.fas',
        blastp = '.snakemake/conda/f1cfc911/bin/blastp'
    conda: ENV
    script: 'scripts/blast_marker_gene.R'

rule novel_marker_genes:
    output: 
        csv = "outputs/GhostKOALA/marker-genes/novel_crumb_marks.csv" 
    input:
        crumb_kegg = 'outputs/GhostKOALA/user_ko_definition.txt',
        bin_kegg = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt'
    conda: ENV
    script: 'scripts/crumb_novel_marker_genes.R'

rule nitro_genes:
    output: csv = 'outputs/GhostKOALA/other/crumb_nitro_genes.csv'
    input:
        info = 'inputs/hu_info.csv',
        crumb_kegg = 'outputs/GhostKOALA/user_ko_definition.txt',
        bin_kegg = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt'
    conda: ENV
    script: 'scripts/crumb_nitrogenase.R'
 
rule plot_fig5a:
    output: 
        pdf = 'outputs/figures/fig5a.pdf',
        png = 'outputs/figures/fig5a.png'
    input:
        info = 'inputs/hu_info.csv',
        crumb_kegg = 'outputs/GhostKOALA/user_ko_definition.txt',
        bin_kegg = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt',
    conda: ENV
    script: 'scripts/fig5a.R'

rule plot_fig5b:
    output:
        pdf = 'outputs/figures/fig5b.pdf',
        png = 'outputs/figures/fig5b.png' 
    input:
        info = "inputs/hu_info.csv",
        crumb_kegg = 'outputs/GhostKOALA/user_ko_definition.txt'
    conda: ENV
    script: 'scripts/fig5b.R'

