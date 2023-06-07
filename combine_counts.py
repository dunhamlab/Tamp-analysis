import sys
# import datetime

# change working directory to where counts files are
workdir = sys.argv[2] + '/merge_files/'
scriptdir = sys.argv[2] + '/Scripts/'

# run w/ python combine_counts.py MN1D_up.txt (from Scripts directory)

samp_file = open(scriptdir + sys.argv[1], 'r')

# first line is generations (header of outfile)
gens = samp_file.readline().strip()

samp_name = sys.argv[1].split('.')[0]

samp_list = []		# list of sample numbers (F1, F4, etc.)
for line in samp_file:
	newsamp = line.strip()
	samp_list.append(newsamp)

file_list = samp_list

line_lists = ['']*len(samp_list)	# list of lineLists from timepoint files
dict_list = [{}]*len(samp_list)		# list of dicts for timepoint files

# date = datetime.datetime.today().strftime('%Y%m%d_')[2:]
outpath = workdir + samp_name + '_counts.txt'
outfile = open(outpath, 'w')

outfile.write(gens + '\n')

num_lines = sum(1 for line in open(workdir + samp_list[0] + '_counts.txt')) + 1

check_list = ['']*len(samp_list)	# list to check genes are in the same order across timepoints
i=0
for item in samp_list:
	filename = workdir + item + '_counts.txt'
	file_list[i] = open(filename, 'r')
	line_lists[i] = file_list[i].readline().strip().split(',')
	check_list[i] =  line_lists[i][0]
	i += 1
	
first_gene = check_list[0]
k=0
done = False
while k < num_lines:
	if all(x==check_list[0] for x in check_list) == True:
		gene_dict = {}		#for summing barcode counts and filtering low count barcodes, sums barcode:count pairs
		n = 0
		for timeList in line_lists:  # for timepoints, sum barcode counts across times in gene_dict
			temp_dict = {}
			gene = timeList[0]
			l = len(timeList)-1
			j = 1
			if l > 1:
				while j < l:
					tag = timeList[j]
					count = int(timeList[j+1])
					temp_dict[tag] = count
					if tag not in gene_dict.keys():
						gene_dict[tag] = count
					elif tag in gene_dict.keys():
						gene_dict[tag] = gene_dict[tag] + count
					j += 2
			dict_list[n] = temp_dict
			line_lists[n] = ''
			n += 1
	
		gene_list = gene_dict.items() # barcode:count sum pairs
		sorted_gene_list = sorted(gene_list, key = lambda tup: tup[1], reverse = True) # sort by count largest to smallest
		print 'length of sorted list = ' + str(len(sorted_gene_list))
		tag_list = []
		for a in sorted_gene_list:
			if a[1] > 5*len(samp_list):		## cutoff is 5 * number of timepoints
				if len(tag_list) < 100:			# add barcodes to tag_list to keep for adding to file
					tag_list.append(a[0])
		print 'number of barcodes = ' + str(len(tag_list))						# number of barcodes for gene
		m = 0
		for barcode in tag_list:
			gene_output = ''
			for Dict in dict_list:
				if barcode in Dict.keys():
					gene_output = gene_output + str(Dict[barcode]) + ','
				else:
					gene_output = gene_output + str(0) + ','
			gene_output_list = gene_output.split(',')
			numlist = []
			for num in gene_output_list[1:(len(gene_output_list)-1)]:
				numlist.append(int(num))
			sums = sum(numlist)
			if sums > 0:
				outfile.write(gene + '_' + str(barcode) + ',' + str(gene_output) + '\n')
			m += 1
		i = 0
		for next in samp_list:
			line_lists[i] = file_list[i].readline().strip().split(',')
			check_list[i] =  line_lists[i][0]
			i += 1
		if check_list[0] == first_gene:
			print 'NOT LOOPING'
			done = True
			break
	else:
		print 'ERROR = wrong order'
		print check_list		
		break
	if done:
		print 'ERROR - not looping'
		break
	k+=1
	
outfile.close()
for item in samp_list:
	item.close()



