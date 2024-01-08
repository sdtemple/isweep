# This script writes SLiM demography files for simulations

import sys
from random import randint
import numpy as np
import pandas as pd

macro, micro, Ne, n, timeSplit, V, L, m, q, rho, gcProp, gcMeanLength, a, b, equalSizes = sys.argv[1:]

sims = pd.read_csv(micro, sep='\t', header=0)
J = sims.shape[0]
sims['FOLDER'] = [macro + '/' + sims.loc[j].MICROEXP + '/' + str(sims.loc[j].SIMNAME) for j in range(J)]

# common setup

# subpopulations, migration
m = int(float(m))
assert m >= 1, 'm is an integer greater than or equal to 1'
q = float(q)
if m > 1:
    qm = str(q / (m - 1))
else:
    qm = str(0.0)
if m > 2:
    print('Warning: are you sure you want to study more than 2 subpopulations?')
rng = [i * 2 for i in range(m, 0, -1)]
unequalSizes = np.random.dirichlet(rng, 1)
us = 'c('
for i in unequalSizes[0]:
    us += str(i)
    us += ','
us = us[:-1]
us += ')'

# recombination, gene conversion
timeSplit = int(float(timeSplit))
RHO = float(rho)
assert RHO >= 0, 'rho is nonnegative'
GCP = float(gcProp)
assert GCP >= 0, 'gcProp is nonnegative'
assert GCP < 1, 'gcProp is not 1'
assert float(gcMeanLength) > 0, 'gcMeanLength is positive'
MODRHO = RHO / (1 - GCP)

# genomic region, position, adaptive variants
L = float(L)
L = L / RHO / 100
L = int(L) - 1
assert L > 0, 'L is greater than zero'
pos = int(L / 2)
L = str(L); pos = str(pos)
V = int(float(V))
assert V >= 1, 'V is integer greater than or equal to 1'
if V > 1:
    print('Warning: are you sure to want to start forward simulations with ' + str(V) + ' adaptive variants?')

# frequency bounds
a=float(a)
b=float(b)
assert (q <= 1.0) and (q >= 0.0), 'q is a proportion'
assert (b <= 1.0) and (b >= 0.0), 'b is a probability'
assert (a <= 1.0) and (a >= 0.0), 'a is a probability'
assert b > a, 'b is greater than a'
a=str(a)
b=str(b)

for j in range(J):
    row = sims.loc[j,]
    s = row.SELCOEF
    tau = row.TIMEMUT
    prefix = str(row.FOLDER) + '/slimulation'

    # selection, age, recombination, gene conversion
    s = str(float(s) * 2)
    tau = int(float(tau))

    # determine variants
    nG = sum([1 for _ in open(Ne)]) - 2
    assert tau > 0, 'tau is positive'
    assert tau < nG, 'tau is less than number of generations in population history file'
    if nG - tau > 1:
        assert V < 2, 'V cannot be positive integer greater than 1 if age of adaptive allele is not end of population history'
    assert timeSplit > 0, 'timeSplit is positive'
    if timeSplit >= nG:
        startSplit = 0
    else:
        startSplit = max(nG - timeSplit,2)
    assert timeSplit > tau, 'age of population substructure is greater than age of adaptive variant'

    # these are for loops in slim script
    # early() and late() statements
    tauf=str(nG-tau)
    taufplus=str(nG-tau+1)
    nG = str(nG)

    # strings for writing files
    if m > 1:
        pops = 'c('
        for i in range(1, m+1):
            if i != m:
                pops += 'p' + str(i) + ', '
            else:
                pops += 'p' + str(i)
        pops += ')'
    else:
        pops = 'p1'
    freqs = 'paste(asString(sim.mutationFrequencies(' + pops + '))'
    for i in range(1, m+1):
        freqs += ', asString(sim.mutationFrequencies(p' + str(i) + '))'
    freqs += ', sep = "\\t")'
    freqrow = 'writeFile("' + prefix + '.freq", ' + freqs + ');\n'
    freqrowa = 'writeFile("' + prefix + '.freq", ' + freqs + ', append = T);\n'

    ##### metaprogramming for SLiM script #####

    f = open(prefix + '.slim', 'w')

    ### initialize ###
    f.write('// set up a simple adaptive simulation\n')
    f.write('initialize() {\n')
    f.write('\tinitializeTreeSeq();\n')
    f.write('\t// fixed adaptive mutation\n')
    f.write('\tinitializeMutationRate(0);\n')
    f.write('\tinitializeMutationType("m1", 0.5, "f", ' + s + ');\n')
    f.write('\n')
    f.write('\t// chromosome\n')
    f.write('\tinitializeGenomicElementType("g1", m1, 1.0);\n')
    f.write('\tinitializeGenomicElement(g1, 0, ' + L + ');\n')
    f.write('\tinitializeRecombinationRate(' + str(MODRHO) + ');\n')
    f.write('\tinitializeGeneConversion(' + gcProp + ', ' + gcMeanLength + ', 1.0);\n')
    f.write('\n')
    f.write('\t// bring in effective population sizes\n')
    f.write('\tdefineGlobal("effSize", asInteger(readFile("' + Ne + '")));\n')
    f.write('\tdefineGlobal("numGen", length(effSize));\n')
    f.write('\tdefineGlobal("itrGen", 0);\n')
    f.write('\n')
    f.write('\t// bring in variant\n')
    f.write('\tdefineGlobal("posVar", ' + pos + ');\n')
    f.write('}\n')
    f.write('\n')

    ### create population ###
    if startSplit < 1:
        f.write('// create population\n')
        f.write('1 {\n')
        f.write('\tdefineConstant("simID", getSeed());\n')
        f.write('\tsubpopCount=' + str(m) + ';\n')
        f.write('\tunequalSizes=' + us + ';\n')
        if m > 1:
            f.write('\tfor (i in 1:subpopCount){\n')
            if equalSizes == '1':
                f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * ' + str(1/m) + ');\n')
                f.write('\t\tsim.addSubpop(i,subpopSize);\n')
            else:
                f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * unequalSizes[i-1]);\n')
                f.write('\t\tsim.addSubpop(i,subpopSize);\n')
            f.write('\t}\n')
        else:
            f.write('\tsubpopSize=asInteger(effSize[itrGen]);\n')
            f.write('\tsim.addSubpop(1,subpopSize);\n')
        if m > 1:
            f.write('\tfor (i in 1:subpopCount){\n')
            f.write('\t\tfor (j in 1:subpopCount){\n')
            f.write('\t\t\tif (i!=j){\n')
            f.write('\t\t\t\tsim.subpopulations[i-1].setMigrationRates(sim.subpopulations[j-1],' + qm + ');\n')
            f.write('\t\t\t}')
            f.write('\t\t}\n')
            f.write('\t}\n')
        else:
            pass
        f.write('\titrGen = itrGen + 1;\n')
        f.write('}\n')
        f.write('\n')
    else:
        f.write('// create population\n')
        f.write('1 {\n')
        f.write('\tdefineConstant("simID", getSeed());\n')
        f.write('\tsubpopSize=asInteger(effSize[itrGen]);\n')
        f.write('\tsim.addSubpop(1,subpopSize);\n')
        f.write('\titrGen = itrGen + 1;\n')
        f.write('}\n')
        f.write('\n')

    ### size updates ###
    if startSplit < 1: # substructure split from the beginning
        f.write('// size updates\n')
        f.write('2:' + nG + ' early() {\n')
        f.write('\tsubpopCount=' + str(m) + ';\n')
        f.write('\tunequalSizes=' + us + ';\n')
        f.write('\tfor (i in 1:subpopCount){\n')
        if equalSizes == '1':
            f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * ' + str(1/m) + ');\n')
        else:
            f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * unequalSizes[i-1]);\n')
        f.write('\t\tsim.subpopulations[i-1].setSubpopulationSize(subpopSize);\n')
        f.write('\t}\n')
        f.write('\titrGen = itrGen + 1;\n')
        f.write('}\n')
        f.write('\n')
    else: # substructure split midway
        # before substructure split
        f.write('// size updates\n')
        f.write('2:' + str(startSplit) + ' early() {\n')
        f.write('\tsubpopSize=' + 'asInteger(effSize[itrGen]);\n')
        f.write('\tsim.subpopulations[0].setSubpopulationSize(subpopSize);\n')
        f.write('\titrGen = itrGen + 1;\n')
        f.write('}\n')
        # substructure split
        f.write(str(startSplit+1) + ' early() {\n')
        f.write('\tsubpopCount=' + str(m) + ';\n')
        if m > 1:
            f.write('\tunequalSizes=' + us + ';\n')
            f.write('\tfor (i in 2:subpopCount){\n')
            if equalSizes == '1':
                f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * ' + str(1/m) + ');\n')
                f.write('\t\tsim.addSubpopSplit(i,subpopSize,1);\n')
            else:
                f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * unequalSizes[i-1]);\n')
                f.write('\t\tsim.addSubpopSplit(i,subpopSize,1);\n')
            f.write('\t}\n')
            if equalSizes == '1':
                f.write('\tsubpopSize=asInteger(effSize[itrGen] * ' + str(1/m) + ');\n')
            else:
                f.write('\tsubpopSize=asInteger(effSize[itrGen] * unequalSizes[0]);\n')
            f.write('\tsim.subpopulations[0].setSubpopulationSize(subpopSize);\n')
            f.write('\tfor (i in 1:subpopCount){\n')
            f.write('\t\tfor (j in 1:subpopCount){\n')
            f.write('\t\t\tif (i!=j){\n')
            f.write('\t\t\t\tsim.subpopulations[i-1].setMigrationRates(sim.subpopulations[j-1],' + qm + ');\n')
            f.write('\t\t\t}\n')
            f.write('\t\t}\n')
            f.write('\t}\n')
        else:
            f.write('\tsubpopSize=asInteger(effSize[itrGen]);\n')
            f.write('\tsim.subpopulations[0].setSubpopulationSize(subpopSize);\n')
        f.write('\titrGen = itrGen + 1;\n')
        f.write('}\n')
        # after substructure split
        f.write(str(startSplit+2) + ':' + nG + ' early() {\n') # ideally startSplit+2 <= nG
        f.write('\tsubpopCount=' + str(m) + ';\n')
        f.write('\tunequalSizes=' + us + ';\n')
        f.write('\tfor (i in 1:subpopCount){\n')
        if equalSizes == '1':
            f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * ' + str(1/m) + ');\n')
        else:
            f.write('\t\tsubpopSize=asInteger(effSize[itrGen] * unequalSizes[i-1]);\n')
        f.write('\t\tsim.subpopulations[i-1].setSubpopulationSize(subpopSize);\n')
        f.write('\t}\n')
        f.write('\titrGen = itrGen + 1;\n')
        f.write('}\n')
        f.write('\n')

    ### draw adaptive mutations ###
    f.write('// draw adaptive mutations\n')
    f.write(tauf + ' late() {\n')
    f.write('\ttarget = sample(' + pops + '.genomes, ' + str(V) + ');\n')
    f.write('\ttarget.addNewDrawnMutation(m1, posVar);\n')
    f.write('\t// save simulation state\n')
    f.write('\tsim.treeSeqOutput("/tmp/slim_" + simID + ".trees");\n')
    f.write('\t// record allele frequency\n')
    f.write('\t' + freqrow)
    f.write('}\n')
    f.write('\n')

    ### make conditional simulation ###
    f.write('// make conditional\n')
    f.write(taufplus + ':' + nG + ' late() {\n')
    f.write('\tcatn(sim.generation);\n')
    f.write('\t' + freqrowa)
    f.write('\tif (sim.countOfMutationsOfType(m1) == 0)\n')
    f.write('\t{\n')
    f.write('\t\tfixed = (sum(sim.substitutions.mutationType == m1) == 1);\n')
    f.write('\t\tif (fixed)\n')
    f.write('\t\t{\n')
    f.write('\t\t\tcatn(simID + ": FIXED - RESTARTING");\n')
    f.write('\t\t\tdeleteFile("' +  prefix + '.freq");\n')
    f.write('\t\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");\n')
    f.write('\t\t\t' + freqrow)
    f.write('\t\t\tsetSeed(rdunif(1, 0, asInteger(2 ^ 62) - 1));\n')
    f.write('\t\t\tdefineGlobal("itrGen",' + tauf + ');\n')
    f.write('\t\t}\n')
    f.write('\t\telse\n')
    f.write('\t\t{\n')
    f.write('\t\t\tcatn(simID + ": LOST - RESTARTING");\n')
    f.write('\t\t\tdeleteFile("' +  prefix + '.freq");\n')
    f.write('\t\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");\n')
    f.write('\t\t\t' + freqrow)
    f.write('\t\t\tsetSeed(rdunif(1, 0, asInteger(2 ^ 62) - 1));\n')
    f.write('\t\t\tdefineGlobal("itrGen",' + tauf + ');\n')
    f.write('\t\t}\n')
    f.write('\t}\n')
    f.write('}\n')
    f.write('\n')

    ### sample current ###
    f.write('// sample current \n')
    f.write(str(int(nG) + 1) + ' early() {\n')
    # draw proportionally from all subpops
    f.write('\tsubpopCount=' + str(m) + ';\n')
    f.write('\tunequalSizes=' + us + ';\n')
    f.write('\tfor (i in 1:subpopCount){\n')
    if equalSizes == '1':
        f.write('\t\tsubpopSize=asInteger(' + str(n) + ' * ' + str(1/m) + ');\n')
        f.write('\t\tsim.subpopulations[i-1].setSubpopulationSize(subpopSize);\n')
    else:
        f.write('\t\tsubpopSize=asInteger(' + str(n) + ' * ' + ' unequalSizes[i-1]);\n')
        f.write('\t\tsim.subpopulations[i-1].setSubpopulationSize(subpopSize);\n')
    f.write('\t}\n')
    f.write('\tsim.chromosome.setRecombinationRate(0);\n')
    f.write('}\n')
    f.write('\n')

    ### output TreeSeq ###
    f.write('// output TreeSeq \n')
    f.write(str(int(nG) + 1) + ' late() {\n')
    f.write('\tcatn(sim.mutationFrequencies(' + pops + '));\n')
    f.write('\tif ((sim.mutationFrequencies(' + pops + ') > ' + a + ') & (sim.mutationFrequencies(' + pops + ') < ' + b + '))\n')
    f.write('\t{\n')
    f.write('\t\t' + freqrowa)
    f.write('\t\tsim.treeSeqOutput("' + prefix + '.trees");\n')
    f.write('\t\tsim.simulationFinished();\n')
    f.write('\t}\n')
    f.write('\telse\n')
    f.write('\t{\n')
    f.write('\t\tcatn(simID + ": OUT OF FREQUENCY BOUNDS - RESTARTING");\n')
    f.write('\t\tdeleteFile("' +  prefix + '.freq");\n')
    f.write('\t\tsim.readFromPopulationFile("/tmp/slim_" + simID + ".trees");\n')
    f.write('\t\t' + freqrow)
    f.write('\t\t\tsetSeed(rdunif(1, 0, asInteger(2 ^ 62) - 1));\n')
    f.write('\t\tdefineGlobal("itrGen",' + tauf + ');\n')
    f.write('\t\tsim.chromosome.setRecombinationRate(' + str(MODRHO) + ');\n')
    f.write('\t}\n')
    f.write('}\n')
    f.close()
