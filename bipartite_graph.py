
import sys
import os
import pickle
import math

starting_directory = str(os.getcwd())
script_path = str(os.path.dirname(os.path.realpath(__file__)))

# Create dictionary for zscores
zscore_infile = open(sys.argv[2], 'r')
zscore_dictionary = {}
for index in zscore_infile:
	index_split = index.split()
	zscore_dictionary[index_split[0]] = index_split[1]

# Read in pickled KO to reaction dictionary
ko_reactionpkl_path = script_path + '/support/ko_reaction.pkl'
ko_dict = pickle.load(open(ko_reactionpkl_path, 'rb'))

# Read in pickled reaction to reaction_mapformula dictionary
reaction_mapformulapkl_path = script_path + '/support/reaction_mapformula.pkl'
reaction_dict = pickle.load(open(reaction_mapformulapkl_path, 'rb'))

infile = open(sys.argv[1], 'r')

infile_name = str(sys.argv[1]).split('/')[-1]
infile_name = infile_name.split('.')[0]

directory = str(os.getcwd()) + '/' + infile_name + '.bipartite.files'
if not os.path.exists(directory):	
	os.makedirs(directory)
os.chdir(directory)
print('Output files located in: ' + directory)

outfile_name = infile_name + '.bipartite.graph'
outfile = open(outfile_name,'w')

# Create file for reporting key errors
errorfile_name = infile_name + '.key_error_log.txt'
errorfile = open(errorfile_name, 'w')

triedCountKO = 0
excludedCountKO = 0
triedCountReact = 0
excludedCountReact = 0
totalIncludedReact = 0

compound_list = []
network_list = []
zscore_list = []

# Nested loops to finally convert the KO list to a directed graph of input and output compounds	
for line in infile:
	current_ko = str(line.split()[1]).strip('ko:')
	#current_ko = str(line.split()[0]).strip('ko:')
	triedCountKO += 1
			
	try:
		reaction_number = ko_dict[current_ko]
	except KeyError:
		errorString = 'WARNING: ' + str(current_ko) + ' not found in KO-to-Reaction dictionary. Omitting.\n'
		errorfile.write(errorString)
		excludedCountKO += 1
		continue # Go to next iteration since this data is necessary
		
	for index in reaction_number:
		triedCountReact += 1
		try:
			reaction_collection = reaction_dict[index]
		except KeyError:
			errorString = 'WARNING: ' + str(index) + ' not found in Reaction-to-Compound dictionary. Omitting.\n'
			errorfile.write(errorString)
			excludedCountReact += 1
			continue
		
		for x in reaction_collection:
			totalIncludedReact += 1
			# Spit reaction input and output as well as the list of compounds with each
			reaction_info = x.split(':')
			input_compounds = reaction_info[0].split('|')
			output_compounds = reaction_info[2].split('|')
						
			for input_index in input_compounds:
				network_list.append(''.join([str(input_index), '\t', str(current_ko), '\n']))
				compound_list.append(str(input_index))
			
			for output_index in output_compounds:
				network_list.append(''.join([str(current_ko), '\t', str(output_index), '\n']))
				compound_list.append(str(output_index))

errorfile.write(''.join(['KOs successfully translated to Reactions: ', str(triedCountKO - excludedCountKO), '\n']))
errorfile.write(''.join(['KOs unsuccessfully translated to Reactions: ', str(excludedCountKO), '\n']))

errorfile.write(''.join(['Reactions successfully translated to Compounds: ', str(triedCountReact - excludedCountReact), '\n']))
errorfile.write(''.join(['Reactions unsuccessfully translated to Compounds: ', str(excludedCountReact)]))

errorfile.close()

network_list = list(set(network_list))
compound_list = list(set(compound_list))

input_zscore_dict = {}
output_zscore_dict = {}
composite_zscore_dict = {}

for index in network_list:
	outfile.write(index)
	edge_info = index.split()
		
	# Output
	if edge_info[0][0] == 'K':
		if not edge_info[1] in output_zscore_dict.keys():
		
			try:
				temp_zscore = zscore_dictionary[edge_info[1]]
			except KeyError:
				temp_zscore = 0
		
			output_zscore_dict[edge_info[1]] = [temp_zscore]
		else:
			output_zscore_dict[edge_info[1]].append(temp_zscore)
		
		# Composite	1
		if not edge_info[1] in composite_zscore_dict.keys():
			
			try:
				temp_zscore = zscore_dictionary[edge_info[1]]
			except KeyError:
				temp_zscore = 0
		
			composite_zscore_dict[edge_info[1]] = [temp_zscore]
		else:
			composite_zscore_dict[edge_info[1]].append(temp_zscore)
		continue
		
		# MIGHT NEED TO REVERSE WHAT I DID HERE
	# Input
	elif edge_info[1][0] == 'K':
		if not edge_info[0] in input_zscore_dict.keys():
			try:
				temp_zscore = zscore_dictionary[edge_info[0]]
			except KeyError:
				temp_zscore = 0

			input_zscore_dict[edge_info[0]] = [temp_zscore]
			
		else:
			input_zscore_dict[edge_info[0]].append(temp_zscore)
		
		# Composite	2
		if not edge_info[0] in composite_zscore_dict.keys():
			try:
				temp_zscore = zscore_dictionary[edge_info[0]]
			except KeyError:
				temp_zscore = 0
				
			composite_zscore_dict[edge_info[0]] = [temp_zscore]
		else:
			composite_zscore_dict[edge_info[0]].append(temp_zscore)

outfile.close()


# create and open composite compound zscore files
inputscorefile_name = infile_name + '.input_zscore.txt'
inputscorefile = open(inputscorefile_name, 'w')

outputscorefile_name = infile_name + '.output_zscore.txt'
outputscorefile = open(outputscorefile_name, 'w')

compositescorefile_name = infile_name + '.composite_zscore.txt'
compositescorefile = open(compositescorefile_name, 'w')

for index in compound_list:

	input_scores = input_zscore_dict[index]
	
	final_score = sum(input_scores) / math.sqrt(len(input_scores))
	inputscorefile.write('\t'.join([index, str(final_score), '\n']))
	
	output_scores = output_zscore_dict[index]
	final_score = sum(output_scores) / math.sqrt(len(output_scores))
	outputscorefile.write('\t'.join([index, str(final_score), '\n']))
	
	composite_scores = composite_zscore_dict[index]
	final_score = sum(composite_scores) / math.sqrt(len(composite_scores))
	compositescorefile.write('\t'.join([index, str(final_score), '\n']))

inputscorefile.close()
outputscorefile.close()
compositescorefile.close()
os.chdir(starting_directory)		
