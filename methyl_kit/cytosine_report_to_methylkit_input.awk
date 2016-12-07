#!/usr/bin/awk

# FUNCTION:
## Turns a cytosine report with the following structure:
## <chr>	<base>	<strand>	<countC>	<countT>	<context>	<exact_context>
##
## to a methylKit input file with the following structure:
## <chrBase>	<chr>	<base>	<strand>	<coverage>	<freqC>	<freqT>
##
## While Filtering every cytosine which do not satisifies the conditions:
## 	- minimum & max coverage

# USAGE EXAMPLE:
## awk -v win_len=100 -v min_cov=2 -v max_cov=100 analysis_type="context_dependent"\ 
## output_files="fileA1_CG.methkit;--file--:fileA1_CHG.methkit;--file--:fileA1_CHH.methkit"\
## -f cytosine_report_to_methylkit_input.awk fileA1.cytosine_report

# ARGS TAKEN:
## win_len 			length of the genomic window
## min_cov			min coverage recquired for a C to be retained within a window
## max_cov			max coverage recquired for a C to be retained within a window
## analysis_type	either 'context_dependent', 'context_independent' or 'both' (determines nb of contexts)
## output_files		string of files joined with ";--infile--:" (a file per context)

# HIERARCHY OF CONTEXTS:
## 1  2   3   4
## CG CHG CHH all_C

BEGIN {

	FS = "\t";

	split(output_files, list_files, ";--outfile--:");

	header = "#<chrBase>\t<chr>\t<base>\t<strand>\t<coverage>\t<freqC>\t<freqT>"

	if ( analysis_type == "context_dependent" ) {

		CG_output = list_files[1];
		CHG_output = list_files[2];
		CHH_output = list_files[3];
		print header > CG_output;
		print header > CHG_output;
		print header > CHH_output;

	} else if ( analysis_type == "context_independent" ) {

		C_output = list_files[1];
		print header > C_output;

	} else {

		CG_output = list_files[1];
		CHG_output = list_files[2];
		CHH_output = list_files[3];
		C_output = list_files[4];
		print header > CG_output;
		print header > CHG_output;
		print header > CHH_output;
		print header > C_output; 

	}

	convert_strand["+"]="F";
	convert_strand["-"]="R";

}

{
	if ( $4 + $5 >= min_cov && $4 + $5 <= max_cov) { 

		chr = $1
		base = $2;
		strand = $3;
		countC = $4;
		countT = $5;
		context = $6;

		coverage = countC + countT;
		freqC = countC / coverage * 100;
		freqT = countT / coverage * 100;

		line = chr "." base "\t" chr "\t" base "\t" convert_strand[strand] "\t" coverage "\t" freqC "\t" freqT;

		if ( analysis_type == "context_dependent" ) {

			if ( context == "CG" ) {
				print line > CG_output;
			} else if ( context == "CHG" ) {
				print line > CHG_output;
			} else {
				print line > CHH_output;
			}

		} else if ( analysis_type == "context_independent" ) {

			print line > C_output;

		} else {

			if ( context == "CG" ) {
				print line > CG_output;
			} else if ( context == "CHG" ) {
				print line > CHG_output;
			} else {
				print line > CHH_output;
			}
			print line > C_output;

		}
	}
}