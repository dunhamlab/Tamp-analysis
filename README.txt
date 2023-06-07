# Tamp-analysis
Analysis scripts for Tamp barseq data

Shell order:
Tamp_analysis
Tamp_combine
Tamp_frequencies
Tamp_plots

Wrapper:
run_qsub.py
	- run with: python run_qsub.py samples.txt analysis.sh
	- will submit jobs qsub analysis.sh lineFromSamplesTxt
	
Shells:
Tamp_analysis.sh
	-UPDATE HEADER AND DIRECTORY
	-run with python run_qsub.py {experiment}_samples.txt Tamp_analysis.sh
	-filename variables might have to be changed depending on the fastq names
	-unzips fastqs, merges read1 and read2, rezips fastqs, removes intermediate .reads files
	-counts barcodes with 190311_count_barcodes.py and 180327_gene_master_list.txt (hard coded location, 
		option for if in script directory)
Tamp_combine.sh
	-UPDATE HEADER AND DIRECTORY
	-run with python run_qsub.py {experiment}_groups.txt Tamp_combine.sh
	-{experiment}_groups.txt file has groups that each have their own corresponding file (below)
	--combine_counts.py
Tamp_frequencies.sh
	-UPDATE HEADER AND DIRECTORY
	-run w/ python run_qsub.py {experiment}_exps.txt Tamp_frequencies.sh
Tamp_plots
	-UPDATE HEADER AND DIRECTORY
	-run w/ qsub Tamp_plots.sh MN_D MN_R DMSOvRadicicol
	-190228_fitness_plots.R
	-190301_plot_chromosomes.R

Scripts:
190311_count_barcodes.py
	-rewrite of AS’s count_barcodes script to be a little faster
		-fixes some dictionary references and reads in barcode dictionary rather than reading and closing file
	-counts barcodes in merge files and creates F*_counts.txt files
	-still has trouble if one gene has a TON of barcodes (like 100K+)
		-in this case I grep searched the merge file to split out that gene from the rest and count that line separately
combine_counts.py
	-combines time points into experiment groups (ex MN1D_up) based on {experiment}_groups.txt
	-filters out barcodes with less that 5*(#timepoints-1), or if t0 count is 0, or if only t0 has reads
	-outfile **_down_counts.txt
190225_calc_frequencies.py
	-adds 1 to all counts and calculates frequencies for every timepoint
	-also stores first line (with generations) in a {exp}_gens.txt file
	-infile = *up/down_counts.txt
	-outfile = *_table.csv and *_gens.txt
190225_calc_log2.py
	-takes log2 ratio of every timepoint to time 0
	-infile = *_table.csv
	-outfile = *_log2.csv
190227_getSlopes.R
	-calls up 190227_CalculateSlopes.R and saves data in {exp}_slopes.csv
190227_CalculateSlopes.R
	-defines function to get slopes for all replicate Tamps
	-chooses either linear regression, or with piecewise linear with one knot at median time based on ANOVA test
	-if piecewise linear, slope of first segment is chosen
	-plots for each gene are put in {exp}_plots directory
	-if #replicates >= 10, fitness is the mean of all reps with standard error
	-if #replicates > 10, fitness is mode of histogram of all reps
190227_remove_NAs.py
	-removes lines with NA (no counts)
190227_dataAnalysis.R
	-makes histogram plots and calculates error cutoff
	-filters out genes with greater than mean+1sd error
	-also normalizes fitness values to the new average of the pool (so new mean is 0)
190228_format_files.py
	-formats files with gene, slope, se, start, stop, chr length, tamp length, etc
	-also uses genes_locations_lengths.tsv
190228_fitness_plots.R
	-various fitness plots/comparisons
190301_plot_chromosomes.R
	-plots for fitness v chromosome coordinate

Support files:
{experiment}_samples.txt
	-list of primer numbers for the experiments
180327_gene_master_list.txt
	-list of all genes possibly in collection
{experiment}_groups.txt
	-file has a list of groups (ex MN1D_up) that each have their own corresponding file (below)
{group}.txt
	-first line is the header for future files with generations
		-ex 'Genes,g0,g2.38,g4.31,g5.61,g6.13,g6.46,g12.77,g19.31,'
	-rest of the lines are the F# to call up F*_counts.txt files
{experiment}_exps.txt
	-the base experiments for combining flask 1&2 and up/down counts
		-ex 'MN_D'
	-should be MN_{exp}

Breakpoint analysis scripts:
tamp-analysis-sep-arms-oct-2018-update-color.R
	-main analysis script for breakpoint analysis
plot-linear-flasso-oct-2018-update-color.R
	-plots data with fused lasso and linear models per chromosome
compare-linear-flasso-sep-arms-oct-2018-update.R
	-Contains functions for variations on fused lasso and linear models
format_breaks.R
	-Determines magnitude of breakpoints identified by fused lasso
	-plots histograms for each experiment
