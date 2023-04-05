import sys
experiments,fileout=sys.argv[1:]
# experiments is *.txt file with filenames of *.tsv files to concat
# fileout is *.tsv file to write to
f=open(fileout,'w')
f.write('MICROEXP\tSIMNAME\tSELCOEF\tTIMEMUT\tSAMPSIZE\n')
with open(experiments,'r') as g:
    for line in g:
        h=open(line.strip(),'r')
        h.readline()
        for line2 in h:
            f.write(line2)
        h.close()
f.close()
