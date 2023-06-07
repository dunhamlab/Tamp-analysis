## called from Tamp_frequencies.sh script
## python 190225_calc_log2.py MN.D /full/dir/Tamp_pool 
import sys
import math
import numpy
from numpy import *

print 'calculating log2'
workdir = sys.argv[2] + '/merge_files/'
scriptdir = sys.argv[2] + '/Scripts/'

expList = sys.argv[1].split('_')

freqFile = open(workdir + sys.argv[1] + '_table.csv', 'r')
gensFile = open(scriptdir + sys.argv[1] + '_gens.txt', 'r')

outFile = open(workdir + sys.argv[1] + '_log2.csv', 'w')

# format generations from gensFile
genList = []
for line in gensFile:
	linelist = line.strip().split(',')
	gens = linelist[2:-1]
	i = 0
	for item in gens:
		number = item[1:]
		gens[i] = str(number)
		i += 1
	genLine = ','.join(gens)
	genList.append(genLine)
	count = len(gens)
gensFile.close()

tlist = ['t']*(count)
glist = ['g']*(count)
numlist = range(1,count+1)
newnum = []
for item in numlist:
	new = str(item)
	newnum.append(new)
timelist = [x+y for x,y in zip(tlist,newnum)]
gslist = [x+y for x,y in zip(glist,newnum)]

firstline = 'gene,mer,replicate,' + ','.join(timelist) + ',' + ','.join(gslist) + '\n'
outFile.write(firstline)


for line in freqFile:
	linelist = line.strip().split(',')
	log2 = []
	for value in linelist[4:]:
		log2.append(float(math.log(float(value)/float(linelist[3]),2)))
	if linelist[2] == 'A':
		Gen = genList[0]
	elif linelist[2] == 'B':
		Gen = genList[1]
	log2str = []
	for val in log2:
		new = str(val)
		log2str.append(new)
	outFile.write(str(linelist[0]) + ',' + str(linelist[1] + ',' + str(linelist[2]) + ',' + ','.join(log2str) + ',' + Gen + '\n'))
	
freqFile.close()
outFile.close()
	
	
	
	
	