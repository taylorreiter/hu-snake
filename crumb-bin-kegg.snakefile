rule download_kegg_uni_crumbs:
    output:        
        ann = 'outputs/hu-crumbs_bin/unitigs/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs_bin/unitigs/GhostKOALA/user.out.top'  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

rule split_kegg_uni_crumbs:
    output:
        ann = dynamic('outputs/hu-crumbs_bin/unitigs/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs_bin/unitigs/GhostKOALA/{crumb_bin}.ko-tax')  
    input: 
        ann = 'outputs/hu-crumbs_bin/unitigs/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs_bin/unitigs/GhostKOALA/user.out.top'  
    run:
        import pandas as pd
        hu = pd.read_table(str(input.ann), header=None)
        hu[0] = hu[0].replace("user:",  "", regex = True)
        hu[6] = hu[0].replace("mhuni_[0-9]{5,6}", "", regex = True)
        hu[6] = hu[6].replace("uni_[0-9]{5,6}", "", regex = True)
        for sample, hu_bin in hu.groupby(6):
            hu_bin.to_csv(f"outputs/hu-crumbs_bin/unitigs/GhostKOALA/{sample}.ko-ann-full.txt", index = False, header = False, sep = "\t")

        hu = pd.read_table(str(input.tax), header=None)
        hu[0] = hu[0].replace("user:",  "", regex = True)
        hu[7] = hu[0].replace("mhuni_[0-9]{5,6}", "", regex = True)
        hu[7] = hu[7].replace("uni_[0-9]{5,6}", "", regex = True)
        for sample, hu_bin in hu.groupby(7):
            hu_bin.to_csv(f"outputs/hu-crumbs_bin/unitigs/GhostKOALA/{sample}.ko-tax", index = False, header = False, sep = "\t")


rule download_kegg_sub_crumbs:
    output: 
        ann = 'outputs/hu-crumbs_bin/subtracts/GhostKOALA/user_ko_defintion.txt',
        tax = 'outputs/hu-crumbs_bin/subtracts/GhostKOALA/user.out.top'  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

rule split_kegg_sub_crumbs:
    output:
        ann = dynamic('outputs/hu-crumbs_bin/subtracts/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs_bin/subtracts/GhostKOALA/{crumb_bin}.ko-tax')  
    input: 
        ann = 'outputs/hu-crumbs_bin/subtracts/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs_bin/subtracts/GhostKOALA/user.out.top'  
    run:
        import pandas as pd
        hu = pd.read_table(str(input.ann), header=None)
        hu[0] = hu[0].replace("user:",  "", regex = True)
        hu[6] = hu[0].replace("mhsub_[0-9]{5,6}", "", regex = True)
        hu[6] = hu[6].replace("sub_[0-9]{5,6}", "", regex = True)
        for sample, hu_bin in hu.groupby(6):
            hu_bin.to_csv(f"outputs/hu-crumbs_bin/subtracts/GhostKOALA/{sample}.ko-ann-full.txt", index = False, header = False, sep = "\t")

        tax = pd.read_table(str(input.tax), header=None)
        tax[0] = tax[0].replace("user:",  "", regex = True)
        tax[7] = tax[0].replace("mhsub_[0-9]{5,6}", "", regex = True)
        tax[7] = tax[7].replace("sub_[0-9]{5,6}", "", regex = True)
        for sample, tax_bin in tax.groupby(7):
            tax_bin.to_csv(f"outputs/hu-crumbs_bin/subtracts/GhostKOALA/{sample}.ko-tax", index = False, header = False, sep = "\t")


rule download_kegg_assembly_crumbs:
    output: 
        ann = 'outputs/hu-crumbs_bin/assembly/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs_bin/assembly/GhostKOALA/user.out.top'  
    shell:'''
    # download kegg annotation files from osf
    # placeholder:
    touch {output.ann}
    touch {output.tax}
    '''

rule split_kegg_assembly_crumbs:
    output:
        ann = dynamic('outputs/hu-crumbs_bin/assembly/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        tax = dynamic('outputs/hu-crumbs_bin/assembly/GhostKOALA/{crumb_bin}.ko-tax')  
    input: 
        ann = 'outputs/hu-crumbs_bin/assembly/GhostKOALA/user_ko_definition.txt',
        tax = 'outputs/hu-crumbs_bin/assembly/GhostKOALA/user.out.top'  
    run:
        import pandas as pd
        hu = pd.read_table(str(input.ann), header=None)
        hu[0] = hu[0].replace("user:",  "", regex = True)
        hu[6] = hu[0].replace("ass_[0-9]{5,6}", "", regex = True)
        for sample, hu_bin in hu.groupby(6):
            hu_bin.to_csv(f"outputs/hu-crumbs_bin/assembly/GhostKOALA/{sample}.ko-ann-full.txt", index = False, header = False, sep = "\t")

        tax = pd.read_table(str(input.tax), header=None)
        tax[0] = tax[0].replace("user:",  "", regex = True)
        tax[7] = tax[0].replace("ass_[0-9]{5,6}", "", regex = True)
        for sample, tax_bin in tax.groupby(7):
            tax_bin.to_csv(f"outputs/hu-crumbs_bin/assembly/GhostKOALA/{sample}.ko-tax", index = False, header = False, sep = "\t")

# UNITIGS ------------------------------------------------

rule plot_kegg_taxonomy_unitigs:
    output: 'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/{crumb_bin}.pdf'
    input: 'outputs/hu-crumbs_bin/unitigs/GhostKOALA/{crumb_bin}.ko-tax'
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
        plot1='outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        plot2b= 'outputs/hu-crumbs_bin/unitigs/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs_bin/unitigs/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs_bin/unitigs/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_unitigs:
    output:
        tsv='outputs/hu-crumbs_bin/unitigs/an-et-al/an_uni.tsv',
        plot='outputs/hu-crumbs_bin/unitigs/an-et-al/an_uni.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs_bin/unitigs/GhostKOALA/{crumb_bin}.ko-ann-full.txt')
    params: 
        title = "Unitigs",
        dir = 'outputs/hu-crumbs_bin/unitigs/GhostKOALA/'
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''

# SUBTRACTS -------------------------------------------------

rule plot_kegg_taxonomy_subtracts:
    output: 'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/{crumb_bin}.pdf'
    input: 'outputs/hu-crumbs_bin/subtracts/GhostKOALA/{crumb_bin}.ko-tax'
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
        plot1='outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        plot2b= 'outputs/hu-crumbs_bin/subtracts/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs_bin/subtracts/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs_bin/subtracts/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_subtracts:
    output:
        tsv='outputs/hu-crumbs_bin/subtracts/an-et-al/an_sub.tsv',
        plot='outputs/hu-crumbs_bin/subtracts/an-et-al/an_sub.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs_bin/subtracts/GhostKOALA/{crumb_bin}.ko-ann-full.txt')
    params: 
        title = "Subtracts",
        dir = 'outputs/hu-crumbs_bin/subtracts/GhostKOALA/'
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''
            
# ASSEMBLIES ------------------------------------------------

rule plot_kegg_taxonomy_assemblies:
    output: 'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/{crumb_bin}.pdf'
    input: 'outputs/hu-crumbs_bin/assembly/GhostKOALA/{crumb_bin}.ko-tax'
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
        plot1='outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-each.pdf',
        plot2a= 'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-all-col-crumb-bin.pdf',
        plot2b= 'outputs/hu-crumbs_bin/assembly/GhostKOALA-plots/kegg-all-col-kingdom.pdf'
    input:
        files = dynamic('outputs/hu-crumbs_bin/assembly/GhostKOALA/{crumb_bin}.ko-ann-full.txt'),
        info = 'sandbox/hu_info.csv'
    params: dir = "outputs/hu-crumbs_bin/assembly/GhostKOALA"
    conda: 'env-kegg.yml'
    shell:'''
    Rscript --vanilla scripts/plot_kegg_anno.R {input.info} {params.dir} {output.plot1} {output.plot2a} {output.plot2b}
    '''

rule plot_an_annotations_assembly:
    output:
        tsv='outputs/hu-crumbs_bin/assembly/an-et-al/an_ass.tsv',
        plot='outputs/hu-crumbs_bin/assembly/an-et-al/an_ass.pdf'
    input:
        an ='sandbox/an-et-al/an_et_al_genes_parsed.tsv',
        files = dynamic('outputs/hu-crumbs_bin/assembly/GhostKOALA/{crumb_bin}.ko-ann-full.txt')
    params: 
        title = "Assembly",
        dir = 'outputs/hu-crumbs_bin/assembly/GhostKOALA/'
    conda: 'env-skimr.yml'
    shell:'''
    Rscript --vanilla scripts/plot_an_et_al.R {input.an} {params.dir} {output.tsv} {params.title} {output.plot}
    '''
