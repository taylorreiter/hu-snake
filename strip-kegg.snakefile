# CRUMB ASSEMBLIES #################################################

rule download_kegg_assembly_crumbs:
    output: 
        ann = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user.out.top'  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

rule split_kegg_assembly_crumbs:
    output:
        ann = dynamic('outputs/hu-crumbs-bin/assembly/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs-bin/assembly/GhostKOALA/{crumb_bin}.ko-tax')  
    input: 
        ann = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user.out.top'  
    run:
        import pandas as pd
        hu = pd.read_table(str(input.ann), header=None)
        hu[0] = hu[0].replace("user:",  "", regex = True)
        hu[6] = hu[0].replace("ass_[0-9]{5,6}", "", regex = True)
        for sample, hu_bin in hu.groupby(6):
            hu_bin.to_csv(f"outputs/hu-crumbs-bin/assembly/GhostKOALA/{sample}.ko-ann-full.txt", index = False, header = False, sep = "\t")

        tax = pd.read_table(str(input.tax), header=None)
        tax[0] = tax[0].replace("user:",  "", regex = True)
        tax[7] = tax[0].replace("ass_[0-9]{5,6}", "", regex = True)
        for sample, tax_bin in tax.groupby(7):
            tax_bin.to_csv(f"outputs/hu-crumbs-bin/assembly/GhostKOALA/{sample}.ko-tax", index = False, header = False, sep = "\t")

# PLOTS

rule plot_kegg_taxonomy_assemblies:
    output: 'outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/{crumb_bin}.pdf'
    input: 'outputs/hu-crumbs-bin/assembly/GhostKOALA/{crumb_bin}.ko-tax'
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

rule plot_kegg_annotations_assembly:
    output: 
        plot1='outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        plot2b= 'outputs/hu-crumbs-bin/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs-bin/assembly/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs-bin/assembly/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_assembly:
    output:
        tsv='outputs/hu-crumbs-bin/assembly/an-et-al/an_ass.tsv',
        plot='outputs/hu-crumbs-bin/assembly/an-et-al/an_ass.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs-bin/assembly/GhostKOALA/{crumb_bin}.ko-ann-full.txt')
    params: 
        title = "Assembly",
        dir = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/'
    conda: 'env-skimr.yml'
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
        blast_nuc = 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc.csv',
        blast_nuc_full = 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_nuc_full.csv',
        blast_aas = 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_aa.csv',
        blast_aas_full = 'outputs/hu-crumbs-bin/assembly/marker-genes/blast_aas_full.csv'
    input:
        info = 'inputs/hu_info.csv',
        crumb_kegg = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt',
        bin_kegg = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt',
        ffn = 'outputs/hu-crumbs-bin/assembly/prokka/all.ffn',
        fna = 'outputs/hu-crumbs-bin/assembly/prokka/all.faa',
        blastn = '.snakemake/conda/b63dce9d/bin/blastn',
        blastp = '.snakemake/conda/b63dce9d/bin/blastp'
    conda: 'env-kegg.yml'
    script: 'scripts/blast_marker_gene.R'

rule plot_fig5a:
    output: 
        pdf = 'outputs/figures/fig5a.pdf',
        png = 'outputs/figures/fig5a.png'
    input:
        info = 'inputs/hu_info.csv',
        crumb_kegg = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt',
        bin_kegg = 'outputs/hu-bins/GhostKOALA/user_ko_definition.txt'
    conda: 'env-kegg.yml'
    script: 'scripts/fig5a.R'

rule plot_fig5b:
    output:
        pdf = 'outputs/figures/fig5b.pdf',
        png = 'outputs/figures/fig5b.png' 
    input:
        crumb_kegg = 'outputs/hu-crumbs-bin/assembly/GhostKOALA/user_ko_definition.txt'
    conda: 'env-kegg.yml'
    script: 'scripts/fig5b.R'
