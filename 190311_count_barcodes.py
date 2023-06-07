#This program will process a file with the yeast barcode in the first column and the bacterial barcode in the second column
#to run: python 140205_count_barcodes.py merge_file 140315_unique_barcodes outfile
import sys
yeast_barcodes = open(sys.argv[2], 'r')		#list of unique barcodes
count_file = open(sys.argv[3], 'w')			#output file

#make a dictionary where the key is the gene and the entry is the deletion collection barcode
tag_dict = {}							#tag:gene
for line in yeast_barcodes:
	gene = line.strip().split()[0]
	tag = line.strip().split()[1]
	tag_dict[tag] = gene

print 'done with dict'

# #make a dictionary where the yeast barcode is the key and the bacterial barcode and the counts per each bacterial barcode is the entry
barcode_dict = {}						#gene:[[bcb,count], [bcb2, count2]]

data_file = open(sys.argv[1], 'r')
data_list = []
i=0
for line in data_file:
	bcy = line.strip().split()[0]
	bcb = line.strip().split()[1]
	new_list = [bcy,bcb]
	data_list.append(new_list)
	if i%100000 == 0:
		print i
	i += 1

for tag in tag_dict:
	tag2_dict = {}
	tag2_string = ''
	for pair in data_list:
		bcy = pair[0]
		bcb = pair[1]
		if tag in bcy:
			if bcb not in tag2_dict:
				tag2_dict[bcb] = 1
			elif bcb in tag2_dict:
				tag2_dict[bcb] = tag2_dict[bcb] +1
	print 'number of bacterial barcodes in ' + str(tag) + '= ' + str(len(tag2_dict))
	
	for item in tag2_dict:
		tag2_string = tag2_string + str(item) + ',' + str(tag2_dict[item]) + ','
	count_file.write(tag_dict[tag] + ',' + tag2_string + '\n')
	
yeast_barcodes.close()
count_file.close()
	
	
	