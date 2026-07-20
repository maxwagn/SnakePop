############################################################

rule popstats:
    input:
        rules.pca.input,
        rules.popgenwindows.input,
        rules.winpca.input,
        rules.manhattan.input,
        rules.heterozygosity.input,
        rules.roh.input,
        rules.snptrees_iqtree.input,
        rules.astral.input,
        rules.dsuite.input,
        rules.twisst.input,
        rules.sfs.input
