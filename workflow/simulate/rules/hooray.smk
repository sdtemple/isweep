wildcard_constraints:
	SIMNAME = '\w+',

rule hooray:
    input:
        fileout='{macro}/{micro}/{seed}/isweep.inference.tsv.gz',
        ranks='{macro}/{micro}/{seed}/isweep.ranks.tsv.gz',
    output:
        hooray='{macro}/{micro}/{seed}/hooray.txt',
        yaml='{macro}/{micro}/{seed}/arguments.yaml',
    params:
        yaml=str(config['CHANGE']['FOLDERS']['YAML']),
    shell:
        """
        touch {output.hooray}
        cp {params.yaml} {output.yaml}
        """
