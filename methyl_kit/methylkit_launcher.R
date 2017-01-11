#!/usr/bin/env Rscript

#############################################################################
# LOAD THE LIBRARIES
#############################################################################
library(optparse);
# Installation
# > install.packages("getopt")
# > library(devtools)
# > install_github("trevorld/optparse", repos=NULL, ref="master", dependencies=TRUE, type="source")
library(plyr); 
library(methylKit); 
library(data.table)
# library(tictoc)


#############################################################################
# READ THE OPTIONS 
#############################################################################
option_list = list(
	make_option("--n_threads", dest="n_threads", action="store"),
	make_option("--name_group_A", dest="name_group_A", action="store"),
	make_option("--name_group_B", dest="name_group_B", action="store"),
	make_option("--files_A", dest="str_files_A", action="store"),
	make_option("--ids_A", dest="str_names_A", action="store"),
	make_option("--files_B", dest="str_files_B", action="store"),
	make_option("--ids_B", dest="str_names_B", action="store"),
	make_option("--context", dest="context", action="store"),
	make_option("--diff", dest="diff", action="store", type="integer"),
	make_option("--qv", dest="qv", action="store", type="double"),
	make_option("--win", dest="win", action="store", type="integer"),
	make_option("--norm", dest="norm", action="store"),
	make_option("--pool", dest="pool", action="store"),
	make_option("--out_dir", dest="out_dir", action="store"),
	make_option("--win_report", dest="win_report", action="store"),
	make_option("--plot", dest="plot", action="store"),
	make_option("--plot_cov", dest="plot_cov", action="store"),
	make_option("--plot_meth", dest="plot_meth", action="store"),
	make_option("--plot_diff", dest="plot_diff", action="store"),
	make_option("--plot_cor", dest="plot_cor", action="store"),
	make_option("--plot_pca", dest="plot_pca", action="store")
);

#############################################################################
# EXTRACT THE OPTIONS 
#############################################################################
opt=parse_args(OptionParser(option_list=option_list));

#############################################################################
# CREATE A LIST OF FILENAMES & IDS (first elements are A, last B)
#############################################################################
files_A=strsplit(opt$str_files_A, split=";--infile--:");
files_B=strsplit(opt$str_files_B, split=";--infile--:");
files.list=as.list(c(files_A[[1]], files_B[[1]]));
names_A=strsplit(opt$str_names_A, split=";--id--:");
names_B=strsplit(opt$str_names_B, split=";--id--:");
ids.list=as.list(c(names_A[[1]], names_B[[1]]));
nA=length(files_A[[1]]);
nB=length(files_B[[1]]);

#############################################################################
# CREATE THE TREATMENT VECTOR
#############################################################################
# This Vector must be of len length(files.list) and assign 0 for each index
# of the list corresponding to the groupA, 1 to the group B
############################################################################# 
treatment=c(rep(0, nA), rep(1, nB));

#############################################################################
# CREATE ONE METHYLRAW OBJECT PER INPUT FILE
#############################################################################
# Example of a methyl row object:
# methylRaw object with 142072 rows
# --------------
#    chr start end strand coverage numCs numTs
# 1 Chr1   101 101      +        1     1     0
# 2 Chr1   102 102      +        1     1     0
# 3 Chr1   107 107      +        1     0     1
# 4 Chr1   108 108      +        2     1     1
# 5 Chr1   109 109      +        3     3     0
# 6 Chr1   110 110      -        1     1     0
# --------------
# sample.id: Weigl_30-119_2_1_val_1.fq_test_bismark_bt2_pe.deduplicated.CX_report 
# assembly: standard 
# context: CG 
# resolution: base 
#############################################################################
meth_objects <- read(files.list, sample.id=ids.list, 
	assembly="", header=FALSE, treatment=treatment, context=opt$context);

#############################################################################
# NORMALIZE COVERAGE 
#############################################################################
# The function normalizes coverage values between samples using a scaling
# factor derived from differences between median of coverage distributions
#############################################################################
if( opt$norm == "True" ) {
	meth_objects=normalizeCoverage(meth_objects, method="median")
}

#############################################################################
# CREATE THE GENOMIC WINDOWS
#############################################################################
# 1. Create the windows for tilling the genome
# 2. Merge all samples to one table by using base-pair locations that are 
# covered in all samples
# 3. Eventually pool the samples of the same group together

#############################################################################
# Example of a meth_table if pool is true: 
#     chr start  end strand coverageA numCsA numTsA coverageB numCsB numTsB
# 1  Chr1   101  200      *      1670   1045    625       110     63     47
# 2  Chr1   501  600      *      1230    122   1108      1022     95    927
# 3  Chr1   801  900      *      1035    584    451       401    228    173
#
# if pool is false, there will be as many coverage & num(C|T)s columns as 
# there are meth_objects
#############################################################################
windows <- tileMethylCounts(meth_objects, win.size=opt$win, step.size=opt$win);
meth_table <- unite(windows, destrand=FALSE);
if( opt$pool == "True" ) {
	meth_table <- pool(meth_table, c(opt$name_group_A, opt$name_group_B));
}

#############################################################################
# RENAME THE COLUMNS BY THE NAME OF THE SAMPLE TO WHICH THEY BELONG
#############################################################################
# if( opt$pool == "True" ) {
# 	names(meth_table)[5:10]=c("coverageA", "numCsA", "numTsA", "coverageB", "numCsB", "numTsB");
# } else {
# 	columns_A=NULL; columns_B=NULL;
# 	for (i in 1:nA) {
# 		columns_A=c(columns_A, c(paste("coverageA", i, sep=""), paste("numCsA", i, sep=""), paste("numTsA", i, sep="")));
# 	}
# 	for (i in 1:nB) {
# 		columns_B=c(columns_B, c(paste("coverageB", i, sep=""), paste("numCsB", i, sep=""), paste("numTsB", i, sep="")));
# 	}
# 	names(meth_table)[5:(5+3*nA-1)]=columns_A;
# 	names(meth_table)[(5+3*nA):(5+3*nA+3*nB-1)]=columns_B;
# }
if( opt$pool == "True" ) {
	names(meth_table)[5:10]=c("tot_cov_A", "cov_Cs_A", "cov_Ts_A", "tot_cov_B", "cov_Cs_B", "cov_Ts_B");
} else {
	columns_A=NULL; columns_B=NULL;
	for (i in 1:nA) {
		columns_A=c(columns_A, c(paste("tot_cov_A", i, sep=""), paste("cov_Cs_A", i, sep=""), paste("cov_Ts_A", i, sep="")));
	}
	for (i in 1:nB) {
		columns_B=c(columns_B, c(paste("tot_cov_B", i, sep=""), paste("cov_Cs_B", i, sep=""), paste("cov_Ts_B", i, sep="")));
	}
	names(meth_table)[5:(5+3*nA-1)]=columns_A;
	names(meth_table)[(5+3*nA):(5+3*nA+3*nB-1)]=columns_B;
}

#############################################################################
# RUN THE DIFFERENTIAL ANALYSIS
#############################################################################
# And further filter the Differentially Methylated Windows based on q-value 
# and percent methylation difference cutoffs
#############################################################################
# Example of a diff_table:
#     chr start  end strand       pvalue      qvalue  meth.diff
# 1  Chr1   101  200      * 0.2660275273 0.281606897 -5.3021230
# 2  Chr1   501  600      * 0.6670567911 0.511504096 -0.6232002
# 3  Chr1   801  900      * 0.9056033384 0.650126808  0.4327346
#############################################################################
diff_table <- calculateDiffMeth(meth_table, num.cores = opt$n_threads);
filtered_diff_table <- diff_table[diff_table$qvalue < opt$qv & abs(diff_table$meth.diff) >= opt$diff,]; 
hypo_diff <- filtered_diff_table[filtered_diff_table$meth.diff < 0,]$meth.diff;
hyper_diff <- filtered_diff_table[filtered_diff_table$meth.diff > 0,]$meth.diff;

############################################################################
# GET THE LIST OF ALL COORDINATES ALLOWED BY THE RETAINED WINDOWS
############################################################################
# win_ranges stores the ranges of coordinates of all retained windows:
# Example:
#       chr start   end
# 2916 Chr4  7201  7300
# 4445 Chr5 62001 62100
#
# win_coord_sq: Convert these ranges to sequences of allowed coordinates 
# win_chr_rep: Without forgetting the corresponding chromosome name 
# win_winstart_rep: Nor the start of the window
# by combining these vectors, we obtain the resulting allowed_coord table:
#
#      chr start   win_key
# 1:  Chr4  7201 Chr4.7201
# 2:  Chr4  7202 Chr4.7201
# ...
# 99: Chr4  7300 Chr4.7201
############################################################################
win_ranges=getData(filtered_diff_table)[,1:3];
win_coord_seq=as.vector(apply(win_ranges, 1, function(x) x[2]:x[3]));
win_chr_rep=as.vector(apply(win_ranges, 1, function(x) rep(x[1], (as.integer(x[3])-as.integer(x[2])+1))));
win_winstart_rep=as.vector(apply(win_ranges, 1, function(x) rep(x[2], (as.integer(x[3])-as.integer(x[2])+1))));
win_key_rep=paste(win_chr_rep, gsub(" ", "", win_winstart_rep), sep=".");
allowed_coord_df=data.frame(win_key=win_key_rep, chr=win_chr_rep, start=win_coord_seq);
allowed_coord=data.table(allowed_coord_df, key=c("chr", "start"));

############################################################################
# DETERMINE THE LIST OF RETAINED WINDOWS KEYS (<chr>.<start>)
############################################################################
win_key_list=paste(win_ranges[,1], win_ranges[,2], sep=".");
nW=length(win_key_list);
win_key_list=data.table(win_key = win_key_list, key="win_key");

############################################################################
# FOR EACH METHOBJECT, RETAIN ONLY CYTOSINES THAT FALL INTO RETAINED WINDOWS
# AND COUNT THEM
############################################################################
# This is achieved by retaining only cytosines whose coordinates 
# (chr, start) match the list of allowed coordinates for each sample
#
# NB : Working with data.tables is the fastest way to proceed in R
# cur_filtered_object example:
#     id_vector  chr start   end strand coverage numCs numTs    win_key
#  1:         1 Chr1 56380 56380      +       26    26     0 Chr1.56301
#
# A second step of the loop consists in counting the nb of retained 
# cytosines within each window. A left outer join is performed between the 
# list of retained windows and the count of cytosines per windows in order 
# to assign a count of "0" for the windows that are potentially empty of
# cytosines in a given sample.

# The process is repeated with an object storing only cytosines present
# in all samples. Common cytosines are retrieved thanks to the function 
# "unite" 
############################################################################
all_cyt_count_table=data.table(matrix(0, nW, (nA+nB+2)));
for(i in 1:(nA+nB)) {
	cur_filtered_object <- data.table(merge(allowed_coord, data.table(getData(meth_objects[[i]]), key=c("chr", "start"))));
	setkey(cur_filtered_object, win_key);
	cur_covered_cyt_count_table <- cur_filtered_object[, .N, keyby="win_key"];
	# Outter join to assign "0" to potentially non covered windows in the given sample i:
	all_cyt_count_table[,i] <- cur_covered_cyt_count_table[win_key_list, on="win_key"]$N;
}
common_cyt_object <- data.table(merge(allowed_coord, data.table(getData(unite(meth_objects)), key=c("chr", "start"))));
common_cyt_count <- common_cyt_object[, .N, by="win_key"];
all_cyt_count_table[,(i+1)] <- common_cyt_count[win_key_list, on="win_key"]$N;
all_cyt_count_table=data.table(cbind(all_cyt_count_table, win_key_list));

#############################################################################
# RENAME THE COLUMNS BY THE NAME OF THE SAMPLE TO WHICH THEY BELONG
#############################################################################
columns_A=NULL; columns_B=NULL;
for (i in 1:nA) {
	columns_A=c(columns_A, paste("nb_Cs_A", i, sep=""))
}
for (i in 1:nB) {
	columns_B=c(columns_B, paste("nb_Cs_B", i, sep=""))
}
names(all_cyt_count_table)[1:nA]=columns_A;
names(all_cyt_count_table)[(nA+1):(nA+nB)]=columns_B;
names(all_cyt_count_table)[(nA+nB+1)]="nb_Cs_inCommon";

#############################################################################
# WRITE THE RESULTING OUTPUT TABLE
#############################################################################
prefix=paste(opt$name_group_A, opt$name_group_B, sep="_vs_");
if( !is.null(opt$win_report) ) {
	merged_table <- data.table(merge(getData(filtered_diff_table), getData(meth_table), 
		by.x=c("chr", "start", "end"), by.y=c("chr", "start", "end"), sort=FALSE));
	merged_table[, c("strand.x", "strand.y") := NULL]; ## drop these two columns
	win_key_list_in_good_order = paste(merged_table$chr, ".", merged_table$start, sep="");
	final_cyt_count_table=all_cyt_count_table[win_key_list_in_good_order, on="win_key"]
	merged_table = data.table(cbind(merged_table, final_cyt_count_table[,1:(nA+nB+1)]))
	fwrite(merged_table, opt$win_report, 
		quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE);
}

#############################################################################
# PLOTS
#############################################################################
if( !is.null(opt$plot) ) {
	pdf(opt$plot, width=12, height=8);

#############################################################################
# Determine the sample code of each meth_objects, (ie sampleA1, A2...)
#############################################################################
	samplecode.list=NULL; color.list=NULL;
	for (i in 1:nA) { 
		samplecode.list[i]=paste("sampleA", i, sep=""); 
		color.list[i]="darkseagreen3";
	}
	for (i in 1:nB) { 
		samplecode.list[nA+i]=paste("sampleB", i, sep=""); 
		color.list[nA+i]="lightsalmon3";
	}

#############################################################################
# For each meth_objects, plot the distribution of coverages
#############################################################################
	if( !is.null(opt$plot_cov) ) {
		par(mfrow=c(max(nA, nB), 2), oma=c(1,1,3,1));
		for (i in 1:length(meth_objects)) {
			cov <- getData(meth_objects[[i]])$coverage;
			hist_info <- hist(cov, plot=FALSE);
			hist_info$counts = round(hist_info$counts/sum(hist_info$counts)*100, digits=1);
			if(i <= nA) { 
				par(mfg=c(i, 1)); 
			} else { 
				par(mfg=c(i-nA, 2)); 
			}
			plot(hist_info, 
				main=samplecode.list[i], 
				cex.main=1, 
				freq=TRUE, 
				col=color.list[i], 
				label=TRUE, 
				xlab="read coverage per retained cytosine",
				ylab="frequency (%)",
				ylim=c(0, max(hist_info$counts))*1.1);
			mtext(ids.list[i], cex=0.6, line=0.5);
		}
		title(paste("Distribution of", opt$context, "coverage"), cex.main=1.3, outer=TRUE, line=1);
	}

#############################################################################
# For each meth_objects, plot the distribution of methylation percentages
#############################################################################
	if( !is.null(opt$plot_meth) ) {
		plot.new();
		par(mfrow=c(max(nA, nB), 2), oma=c(1,1,3,1));
		# For each meth_objects, 
		for (i in 1:length(meth_objects)) {
			cur_obj = getData(meth_objects[[i]])[,6:7];
			meth = cur_obj[,1]/rowSums(cur_obj[,1:2])*100;
			hist_info <- hist(meth, plot=FALSE);
			if(i <= nA) { 
				par(mfg=c(i, 1)); 
			} else { 
				par(mfg=c(i-nA, 2)); 
			}
			#getMethylationStats(meth_objects[[i]], plot=TRUE, both.strands=FALSE);
			plot(hist_info, 
				main=samplecode.list[i], 
				cex.main=1,  
				col=color.list[i], 
				label=FALSE, 
				xlab="methylation percentage per retained cytosine",
				ylab="count",
				ylim=c(0, max(hist_info$counts))*1.1)
			mtext(ids.list[i], cex=0.6, line=0.5);
		}
		title(paste("Distribution of cytosine methylation % in", opt$context, "context"), cex.main=1.3, outer=TRUE, line=1);
	}

#############################################################################
# Plot the distribution of methylation differences of the windows, taking
# group A as the reference
#############################################################################
	if( !is.null(opt$plot_diff) ) {
		par(mfrow=c(1,1), oma=c(3,8,4,8));
		all_hist_info <- hist(diff_table$meth.diff, breaks=seq(-100,100,1), plot=FALSE);
		ymax=max(all_hist_info$count);
		plot(all_hist_info, 
			main=NA,
			col="gray85", 
			xlab="methylation difference (%)",
			ylab="count",
			xlim=c(-100, 100),
			ylim=c(0, ymax*1.1));

		title(paste("Distribution of methylation difference of genomic windows\n", prefix, "in", opt$context, "context"),
			cex.main=1.2, line=2);
		abline(v=c(-opt$diff, opt$diff), lty="dashed", lwd=2);
		text(-100, ymax*1.1, pos=4, labels ="Hypo methylated area", cex=0.9);
		text(0, ymax*1.1, labels ="Filtered area", cex=0.9);
		text(100, ymax*1.1, pos=2, labels ="Hyper methylated area", cex=0.9);
		legend("right", legend=c("all windows", "retained hypo windows", "retained hyper windows"), 
			fill=c("gray85", "lightblue1", "lightcoral"), cex=0.8);
		
		hypo_hist_info <- hist(hypo_diff, plot=FALSE, breaks=seq(-100,100,1));
		lines(hypo_hist_info, 
			main=NA,
			col="lightblue1");

		hyper_hist_info <- hist(hyper_diff, plot=FALSE, breaks=seq(-100,100,1));
		lines(hyper_hist_info, 
			main=NA,
			col="lightcoral");
	}

#############################################################################
# Plot correlation of methylation % distribution between samples
#############################################################################
	if( !is.null(opt$plot_cor) ) {
		getCorrelation(meth_table, plot=TRUE);
	}

#############################################################################
# Plot PCA 
#############################################################################$
	if( !is.null(opt$plot_pca) ) {
		PCASamples(meth_table, comp=c(1,2));
	}
	junk=dev.off();
}