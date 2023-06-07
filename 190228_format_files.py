import sys
# args: group dir

group = sys.argv[1]
workdir = sys.argv[2] + '/merge_files/'
scriptdir = sys.argv[2] + '/Scripts/'

reffile = open(scriptdir + 'genes_locations_lengths.tsv','r')

infile = open(workdir + group + '_slopes_mean_norm.csv', 'r')
outfile = open(workdir + group + '_slopes_se_norm.csv', 'w')

firstline = 'genes,number,chr,' + group + ',' + group + '.se,start,stop,chrlength,arm,Tamp.length'
outfile.write(firstline + '\n')

infile.readline()

alphadict = {'A':1, 'B':2, 'C':3, 'D':4, 'E':5, 'F':6, 'G':7, 'H':8, 'I':9, 'J':10, 'K':11, 'L':12, 'M':13, 'N':14, 'O':15, 'P':16}
chrlendict = {'A':230218, 'B':813184, 'C':316620, 'D':1531933, 'E':576874, 'F':270161, 'G':1090940, 'H':562643, 'I':439888, 'J':745751, 'K':666816, 'L':1078177, 'M':924431, 'N':784333, 'O':1091291, 'P':948066}

# make reference dictionary
refdict = {}
for line in reffile:
	linelist = line.strip().split('\t')
	refdict[linelist[0]] = [linelist[2], linelist[3], linelist[4]]

i = 0
for line in infile:
	list = line.strip().split(',')
	genes = list[0]
	print genes
	if genes in refdict:
		slope = list[1]
		se = list[2]
		chrID = genes[1]
		chr = alphadict[chrID]
		ref = refdict[genes]
		numstart = ref[0]
		stop = ref[1]
		chrlen = chrlendict[chrID]
		arm = genes[2]
		if arm == 'L':
			Tamplen = numstart
		elif arm == 'R':
			Tamplen = int(chrlen) - int(stop)
		newline = str(genes) + ',' + str(numstart) + ',' + str(chr) + ',' + str(slope) + ',' + str(se) + ',' + str(numstart) + ',' + str(stop) + ',' + str(chrlen) + ',' + str(arm) + ',' + str(Tamplen)
		outfile.write(newline + '\n')
	else:
		print 'NOT FOUND'
		i += 1
print str(i) + ' genes not found'		

reffile.close()
infile.close()
outfile.close()