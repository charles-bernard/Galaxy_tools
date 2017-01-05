#!/usr/bin/awk

# FUNCTION:
## Identify differentially methylated regions from a list of differentially methylated rwindows: 
## <chr> <base> <end> <meth.diff> <minPooledCov> <minCytCount> <ratio_commonCytCount_maxCytCount>
##
## dmr file:
## <chr> <start> <end> <meth_diff> <minCytCount> <max_nb_Cs/nb_Cs_inCommon>
##
## While Filtering every window which do not satisifies the conditions:
##  - minimum coverage
##  - min nb of cytosines of a given context
##  - min ratio max nb cytosine / nb common cytosines

# USAGE EXAMPLE:
## awk -v outdir="outdir" -v name="name" -v context="CG" -v min_C=3 -v min_cov=10 -v min_ratio=0.4 -v n_gap=1
## -f join_dmw_into_dmr.awk statistic_file

# ARGS TAKEN:
## win              win_size
## min_C            recquired nb of cytosines
## min_cov			min coverage recquired for a window to be retained
## min_ratio        ratio recquired

# HIERARCHY OF CONTEXTS:
## 1  2   3   4
## CG CHG CHH all_C

BEGIN {

	FS = "\t";

	old_chr["hypo"] = "";
	old_start["hypo"] = -1;
	old_end["hypo"] = -1;
	nW["hypo"] = 1;
	outfile_std["hypo"] = outdir "/hypo_differenceDMRs.txt";
	outfile_fig["hypo"] = outdir "/hypo_differenceDMRs_for_figure_and_annotation.txt";
	outfile_gff["hypo"] = outdir "/hypo_differenceDMRs.gff";

	old_chr["hyper"] = "";
	old_start["hyper"] = -1;
	old_end["hyper"] = -1;
	nW["hyper"] = 1;
	outfile_std["hyper"] = outdir "/hyper_differenceDMRs.txt";
	outfile_fig["hyper"] = outdir "/hyper_differenceDMRs_for_figure_and_annotation.txt";
	outfile_gff["hyper"] = outdir "/hyper_differenceDMRs.gff";

	header_std = "#chr\tstart\tend\tmeth.diff\tminCytCount\tratio_commonCytCount_maxCytCount";
	header_fig = "#chr\tstart\tend";
	header_gff = "#chr\tsource\tfeature\tstart\tend\tmeth.diff\tstrand\tframe\tattribute";
	print header_std > outfile_std["hypo"];
	print header_std > outfile_std["hyper"];
	print header_fig > outfile_fig["hypo"];
	print header_fig > outfile_fig["hyper"];
	print header_gff > outfile_gff["hypo"];
	print header_gff > outfile_gff["hyper"];

	hexadec_color["CG"] = "4C9001A"; #Apple Green
	hexadec_color["CHG"] = "495CFF"; #Blue
	hexadec_color["CHH"] = "1D702D"; #Dark Green
	hexadec_color["C"] = "4c4c4c";   #Grey
}

NR > 1 {

	if (NR==2) { win = $3 - $2 + 1 }

	if($4/100 > 0) { 
		dmr_type = "hyper";
	} else { 
		dmr_type = "hypo";
	}

	meth_diff = $4/100;
	cov = $5;
	cyt_count = $6;
	if($7 == "") { ratio = 0 } else { ratio = $7 }

	if( cov >= min_cov && cyt_count >= min_C && ratio >= min_ratio ) {


		chr = $1;
		start = $2;
		end = $3;
		
		if(chr != old_chr[dmr_type]) {

			if(old_end[dmr_type] > 0) {

				line_std = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t" sum_cyt_count[dmr_type] "\t" sum_ratio[dmr_type]*1/nW[dmr_type];
				line_fig = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type];
				line_gff = old_chr[dmr_type] "\t" name "\t" name "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t.\t.\t" "context=" context ";color=" hexadec_color[context];
				print line_std > outfile_std[dmr_type];
				print line_fig > outfile_fig[dmr_type];
				print line_gff > outfile_gff[dmr_type];

			}
			nW[dmr_type] = 1;
			old_chr[dmr_type] = chr;
			old_start[dmr_type] = start;
			old_end[dmr_type] = end;
			sum_meth_diff[dmr_type] = meth_diff;
			sum_cyt_count[dmr_type] = cyt_count;
			sum_ratio[dmr_type] = ratio;

		} else {

			if(start >= old_end[dmr_type]+1 && start <= old_end[dmr_type] + (n_gap * win) + 1) {

				nW[dmr_type]++;
				old_end[dmr_type] = end;
				sum_meth_diff[dmr_type] += meth_diff;
				sum_cyt_count[dmr_type] += cyt_count;
				sum_ratio[dmr_type] += ratio;

			} else {

				if(old_end[dmr_type] > 0) {

					line_std = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t" sum_cyt_count[dmr_type] "\t" sum_ratio[dmr_type]*1/nW[dmr_type];
					line_fig = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type];
					line_gff = old_chr[dmr_type] "\t" name "\t" name "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t.\t.\t" "context=" context ";color=" hexadec_color[context];
					print line_std > outfile_std[dmr_type];
					print line_fig > outfile_fig[dmr_type];
					print line_gff > outfile_gff[dmr_type];
				
				}
				nW[dmr_type]=1;
				old_start[dmr_type] = start;
				old_end[dmr_type] = end;
				sum_meth_diff[dmr_type] = meth_diff;
				sum_cyt_count[dmr_type] = cyt_count;
				sum_ratio[dmr_type] = ratio;

			}
		}
	}
}

END {

	dmr_type = "hypo";
	line_std = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t" sum_cyt_count[dmr_type] "\t" sum_ratio[dmr_type]*1/nW[dmr_type];
	line_fig = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type];
	line_gff = old_chr[dmr_type] "\t" name "\t" name "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t.\t.\t" "context=" context ";color=" hexadec_color[context];
	print line_std > outfile_std[dmr_type];
	print line_fig > outfile_fig[dmr_type];
	print line_gff > outfile_gff[dmr_type];


	dmr_type = "hyper";
	line_std = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t" sum_cyt_count[dmr_type] "\t" sum_ratio[dmr_type]*1/nW[dmr_type];
	line_fig = old_chr[dmr_type] "\t" old_start[dmr_type] "\t" old_end[dmr_type];
	line_gff = old_chr[dmr_type] "\t" name "\t" name "\t" old_start[dmr_type] "\t" old_end[dmr_type] "\t" sum_meth_diff[dmr_type]*1/nW[dmr_type] "\t.\t.\t" "context=" context ";color=" hexadec_color[context];
	print line_std > outfile_std[dmr_type];
	print line_fig > outfile_fig[dmr_type];
	print line_gff > outfile_gff[dmr_type];

}
