import os
import sys

# python run_qsub.py samplelist.txt analysis.sh
## change directory before running

sampfile = open(sys.argv[1], 'r')
samplist = sampfile.readlines()
Directory = '/net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/Scripts/'

for line in samplist:
	os.system('qsub ' + Directory + sys.argv[2] + ' ' + line.strip())
	
sampfile.close()