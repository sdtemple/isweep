wildcard_constraints:
	SIMNAME = '\w+',

rule record:
	input:
		'{macro}/{micro}/{seed}/scan.chr1.windowed.tsv.gz',
		'{macro}/{micro}/{seed}/short.chr1.ibd.gz',
		'{macro}/{micro}/{seed}/long.chr1.ibd.gz'
	output:
		yaml='{macro}/{micro}/{seed}/arguments.yaml'
	shell:
		'cp {config[CHANGE][FOLDERS][YAML]} {output.yaml}'

rule tszip:
	input:
		trees='{macro}/{micro}/{seed}/slimulation.trees',
		yaml='{macro}/{micro}/{seed}/arguments.yaml',
	output:
		tsz='{macro}/{micro}/{seed}/slimulation.trees.tsz'
	shell:
		'tszip {input.trees}'
