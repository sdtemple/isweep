wildcard_constraints:
	SIMNAME = '\w+',

rule hooray:
    input:
        fileout='{macro}/{micro}/{seed}/isweep.inference.tsv',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
        trank='{macro}/{micro}/{seed}/isweep.rank.true.tsv',
        # clues='{macro}/{micro}/{seed}/clues.post.npy',
    output:
        yaml='{macro}/{micro}/{seed}/arguments.yaml',
    params:
        yaml=str(config['CHANGE']['FOLDERS']['YAML']),
    shell:
        """
        cp {params.yaml} {output.yaml}
        """
