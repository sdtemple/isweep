wildcard_constraints:
	SIMNAME = '\w+',

rule record:
	input:
		'{macro}/{micro}/hapibd.results.tsv',
		'{macro}/{micro}/ibdends.results.tsv',
	output:
		yaml='{macro}/{micro}/arguments.yaml'
	shell:
		'cp ../config/simstudies.yaml {output.yaml}'

rule tszip:
	input:
		trees='{macro}/{micro}/slimulation.trees',
		yaml='{macro}/{micro}/arguments.yaml',
	output:
		tsz='{macro}/{micro}/slimulation.trees.tsz'
	shell:
		'tszip {input.trees}'
