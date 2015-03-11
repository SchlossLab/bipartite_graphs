
import sys
import os

exp_data = open(sys.argv[1], 'r')
probeIDs = open(sys.argv[2], 'r')
koIDs = open(sys.argv[3], 'r')

infile_name = str(sys.argv[1]).split('/')[-1]
infile_name = infile_name.split('.')[0]

outfile_name = infile_name + '.weighted_KOs.tsv'
outfile = open(outfile_name,'w')

probe_dictionary = {}
for index in probeIDs:
	probe_dictionary[index.split()[1]] = str(index.split()[0])
	
ko_dictionary = {}
for index in koIDs:
	index_split = index.split()
	ko_dictionary[index_split[0].strip('cdf:')] = str(index_split[1].strip('ko:'))

probeIDs.close()
koIDs.close()

gene = 'HELLO WORLD'


# Loop to translate probes
for index in exp_data:

	index_split = index.split()
	
	try:
		gene = probe_dictionary[str(index_split[0])]
	except KeyError:
		print('Probe translation error: ' + str(index_split[0]) + ' included as probe')
		ko = str(index_split[0])
		continue
	
	try:
		ko = ko_dictionary[gene]
	except KeyError:
		print('KO translation error: ' + str(gene) + ' included as gene')
		ko = gene
		continue
	
	outfile.write('\t'.join([ko, str(index_split[1])]))
	outfile.write('\n')
	
exp_data.close()
outfile.close()