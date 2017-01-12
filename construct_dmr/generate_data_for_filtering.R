#!/usr/bin/env Rscript

#############################################################################
# LOAD THE LIBRARIES
#############################################################################
library(optparse);
# Installation
# > install.packages("getopt")
# > library(devtools)
# > install_github("trevorld/optparse", repos=NULL, ref="master", dependencies=TRUE, type="source")
library(data.table);
library(matrixStats);

#############################################################################
# READ THE OPTIONS 
#############################################################################
option_list = list(
	make_option("--infile", dest="input_file", action="store"),
	make_option("--outfile", dest="output_file", action="store"),
	make_option("--nA", dest="nA", action="store", type="integer"),
	make_option("--nB", dest="nB", action="store", type="integer"),
	make_option("--pool", dest="pool", action="store")
);

#############################################################################
# EXTRACT THE OPTIONS 
#############################################################################
opt=parse_args(OptionParser(option_list=option_list));

#############################################################################
# READ THE INPUT FILE
#############################################################################
input_table=fread(opt$input_file, header=TRUE, sep="\t");

#############################################################################
# DETERMINE THE INDEXES FOR A|B TOT COVERAGES, POOL IF NECESSARY AND 
# EXTRACT THE MINIMUM POOLED COVERAGE
#############################################################################
if (opt$pool == "True") {
	pooled_cov_A = input_table[,7];
	pooled_cov_B = input_table[,10];
	cov_mat = matrix(c(pooled_cov_A, pooled_cov_B, nrow=length(pooled_cov_A), ncol=2));
} else {
	index_cov_A = NULL;
	index_cov_B = NULL;
	first=7
	for(i in 1:opt$nA) {
		index_cov_A = c(index_cov_A, (first + (i - 1) * 3));
	}
	first=first + i*3;
	for(i in 1:opt$nB) {
		index_cov_B = c(index_cov_B, (first + (i - 1) * 3));
	}
	cov_A_mat = as.matrix(input_table[,index_cov_A, with=FALSE]);
	cov_B_mat = as.matrix(input_table[,index_cov_B, with=FALSE]);
	cov_mat = cbind(cov_A_mat, cov_B_mat);
}
min_cov = rowMins(cov_mat);

#############################################################################
# DETERMINE THE INDEXES FOR THE CYTOSINE COUNTS, COMPUTE THE RATIO 
# (COMMON CYTOSINE COUNT IN THIS ROW) / MAX(CYTOSINE COUNT IN A ROW)
#############################################################################
n_col=ncol(input_table);
index_cyt_count=seq((n_col - (opt$nA + opt$nB)), n_col);
cyt_count_mat = as.matrix(input_table[,index_cyt_count, with=FALSE]);
last = dim(cyt_count_mat)[2];
min_cyt_count = rowMins(cyt_count_mat[,1:(last-1)]);
ratio = cyt_count_mat[,last] / rowMaxs(cyt_count_mat);

#############################################################################
# CREATE THE STATISTICS TABLE
#############################################################################
statistics_table = data.table(input_table[,c(1,2,3,6)], 
	min_cov=min_cov, 
	minCytCount=min_cyt_count, 
	ratio_commonCytCount_maxCytCount=ratio);

fwrite(statistics_table, opt$output_file, 
 	quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
