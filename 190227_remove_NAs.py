import sys
group = sys.argv[1]
path = sys.argv[2]
filename = path + '/' + group + '_slopes.csv'


infile = open(filename , 'r')

outfile = open(str(path + '/' + group + '_slopes_remove_NAs.csv'), 'w')

for line in infile:
	line = line.replace('"','')
	linelist = line.strip().split(',')
	newline = ','.join(linelist)
	if 'NA' in linelist:
		print 'NA'
	elif 'NA' not in linelist:
		print newline
		outfile.write(newline + '\n')
	
infile.close()
outfile.close()
