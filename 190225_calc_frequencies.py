## called from Tamp_frequencies.sh script
## python 190225_calc_frequencies.py MN1D /full/dir/Tamp_pool 
import sys

workdir = sys.argv[2] + '/merge_files/'
scriptdir = sys.argv[2] + '/Scripts/'

expList = sys.argv[1].split('_')

up1filename = workdir + expList[0] + '1' + expList[1] + '_up_counts.txt'
down1filename = workdir + expList[0] + '1' + expList[1] + '_down_counts.txt'
up2filename = workdir + expList[0] + '2' + expList[1] + '_up_counts.txt'
down2filename = workdir + expList[0] + '2' + expList[1] + '_down_counts.txt'

up1file = open(up1filename, 'r')
down1file = open(down1filename, 'r')
up2file = open(up2filename, 'r')
down2file = open(down2filename, 'r')

outfile = open(workdir + sys.argv[1] + '_table.csv', 'w')
# make 
outGens = open(scriptdir + sys.argv[1] + '_gens.txt', 'w')

fileList = [up1file, down1file, up2file, down2file]

line = ''
n=0
for file in fileList:
	line = file.readline()
	if n == 0:
		outGens.write(line)
	elif n == 2:
		outGens.write(line)
	n+=1
outGens.close()
timeNumList = line.strip().split(',')
timeNum = len(timeNumList) - 2

replicates = ['A', 'A', 'B', 'B']
totals = [[0]*timeNum]*4

# get total counts for timepoints
print 'calculating total reads'
i=0
for file in fileList:
	print i
	for line in file:
		linelist = line.strip().split(',')
		if float(linelist[1]) != 0:
			counts = []
			for item in linelist[1:-1]:
				counts.append(float(item)+1)
			totals[i] = [x+y for x,y in zip(totals[i], counts)]
	i+=1
print totals

for file in fileList:
	file.close()

# determine frequencies and print to outfile
print 'calculating frequencies'
up1file = open(up1filename, 'r')
down1file = open(down1filename, 'r')
up2file = open(up2filename, 'r')
down2file = open(down2filename, 'r')
	
fileList = [up1file, down1file, up2file, down2file]

j=0
for file in fileList:
	file.readline()
	for line in file:
		linelist = line.strip().split(',')
		gene = linelist[0].split('_')[0]
		tag = linelist[0].split('_')[-1]
		if float(linelist[1]) != 0:
			counts = []
			for item in linelist[1:-1]:
				counts.append(float(item)+1)
			freq = [x/y for x,y in zip(counts, totals[j])]
			freqLine = ','.join(str(i) for i in freq)
			outfile.write(gene + ',' + tag + ',' + replicates[j] + ',' + freqLine + '\n')
	j+=1

outfile.close()

for file in fileList:
	file.close()








