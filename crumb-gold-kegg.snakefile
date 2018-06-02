rule download_kegg_uni:
    output: 
        ann = dynamic('outputs/hu-crumbs_gold/unitigs/GhostKOALA/{crumbs_gold}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs_gold/unitigs/GhostKOALA/{crumbs_gold}.ko-tax')  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

rule download_kegg_sub:
    output: 
        ann = dynamic('outputs/hu-crumbs_gold/subtracts/GhostKOALA/{crumbs_gold}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs_gold/subtracts/GhostKOALA/{crumbs_gold}.ko-tax')  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''
    
rule download_kegg_assembly:
    output: 
        ann = dynamic('outputs/hu-crumbs_gold/assembly/GhostKOALA/{crumbs_gold}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs_gold/assembly/GhostKOALA/{crumbs_gold}.ko-tax')  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

# UNITIGS ------------------------------------------------

rule plot_kegg_taxonomy_unitigs:
    output: 'outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/{crumbs_gold}.pdf'
    input: 'outputs/hu-crumbs_gold/unitigs/GhostKOALA/{crumbs_gold}.ko-tax'
    run:
        import pandas as pd
        import matplotlib
        import matplotlib.pyplot as plt
        matplotlib.rcParams['figure.figsize'] = [17, 8]
        
        hu = pd.read_table(str(input), header=None)
        hu = hu.loc[hu[6] > 100] # recommended by GhostKOALA
        hu = hu[~hu[4].str.contains("candidate")] # remove things that others were uncertain about
        hu = hu[~hu[4].str.contains("Candidatus")]
        
        tax = hu.groupby(4)
        tax = tax.filter(lambda x: len(x) > (hu.shape[0]/1500)) # dynamically adjust values on x axis
        
        plot = tax[4].value_counts().plot(kind = 'bar')
        plt.tight_layout()
        plot.get_figure().savefig(str(output), format='pdf') 

rule plot_kegg_annotations_unitigs:
    output: 
        plot1='outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/kegg-all-col-crumbs-gold.pdf',
        plot2b= 'outputs/hu-crumbs_gold/unitigs/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs_gold/unitigs/GhostKOALA/{crumbs_gold}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs_gold/unitigs/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_unitigs:
    output:
        tsv='outputs/hu-crumbs_gold/unitigs/an-et-al/an_uni.tsv',
        plot='outputs/hu-crumbs_gold/unitigs/an-et-al/an_uni.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs_gold/unitigs/GhostKOALA/{crumbs_gold}.ko-ann-full.txt')
    params: 
        title = "Unitigs",
        dir = 'outputs/hu-crumbs_gold/unitigs/GhostKOALA/'
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''

# SUBTRACTS -------------------------------------------------

rule plot_kegg_taxonomy_subtracts:
    output: 'outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/{crumbs_gold}.pdf'
    input: 'outputs/hu-crumbs_gold/subtracts/GhostKOALA/{crumbs_gold}.ko-tax'
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

rule plot_kegg_annotations_subtracts:
    output: 
        plot1='outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/kegg-all-col-crumbs-gold.pdf',
        plot2b= 'outputs/hu-crumbs_gold/subtracts/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs_gold/subtracts/GhostKOALA/{crumbs_gold}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs_gold/subtracts/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_subtracts:
    output:
        tsv='outputs/hu-crumbs_gold/subtracts/an-et-al/an_sub.tsv',
        plot='outputs/hu-crumbs_gold/subtracts/an-et-al/an_sub.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs_gold/subtracts/GhostKOALA/{crumbs_gold}.ko-ann-full.txt')
    params: 
        title = "Subtracts",
        dir = 'outputs/hu-crumbs_gold/subtracts/GhostKOALA/'
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''
            
# ASSEMBLIES ------------------------------------------------

rule plot_kegg_taxonomy_assemblies:
    output: 'outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/{crumbs_gold}.pdf'
    input: 'outputs/hu-crumbs_gold/assembly/GhostKOALA/{crumbs_gold}.ko-tax'
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
        plot1='outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/kegg-all-col-crumbs-gold.pdf',
        plot2b= 'outputs/hu-crumbs_gold/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs_gold/assembly/GhostKOALA/{crumbs_gold}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs_gold/assembly/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_assembly:
    output:
        tsv='outputs/hu-crumbs_gold/assembly/an-et-al/an_ass.tsv',
        plot='outputs/hu-crumbs_gold/assembly/an-et-al/an_ass.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs_gold/assembly/GhostKOALA/{crumbs_gold}.ko-ann-full.txt')
    params: 
        title = "Assembly",
        dir = 'outputs/hu-crumbs_gold/assembly/GhostKOALA/'
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''
