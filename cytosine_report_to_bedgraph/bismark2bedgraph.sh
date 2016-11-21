#!/bin/bash

# read the options
TEMP=`getopt -o e::i:co: --long epi::,infile_cov:,context,output_dir:,tool_dir:,tdf,igv_genome: -n 'report' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
	case "$1" in
		-e | --epi )
			case "$2" in
				"" ) epi="current_job"; shift 2 ;;
				*) epi=$2; shift 2 ;;
			esac ;;
		-i | --infile_cov )
			case "$2" in
				*) infile_cov=$2; shift 2 ;;
			esac ;;
		-c | --context ) 
			context=true; shift ;;
		--tdf )
			tdf=true; shift ;;
		--igv_genome )
			case "$2" in
				*) igv_genome=$2; shift 2 ;;
			esac ;;
		-o | --output_dir )
			case "$2" in
				*) output_dir=$2; shift 2 ;;
			esac ;;
		--tool_dir )
			case "$2" in
				*) tool_dir=$2; shift 2 ;;
			esac ;;
		-- ) shift; break ;;
		* ) echo "Internal error!"; exit 1 ;;
	esac
done

# define outputs according to the options
if [[ "$context" = true ]]; then
	n="4"
	output_types=("CG" "CHG" "CHH" "coverage");
	context_list=("CG" "CHG" "CHH" ".*");
	iscoverage=("false" "false" "false" "true");
else
	n="2"
	output_types=("CXX" "coverage");
	context_list=(".*" ".*");
	iscoverage=("false" "true");
fi

for (( i=0; i<$n; i++)); do
	prefix="$output_dir""/""$epi""_""${output_types[$i]}";
	bedgraph_list[$i]="$prefix"".bedgraph";
	if [[ "$tdf" = true ]]; then 
		tdf_list[$i]="$prefix"".tdf"; 
	fi
done

# process function for the conversion from report to bedgraph
awk_process() {
	tool_dir=$1;
	report_infile=$2;
	bedgraph_outfile=$3;
	output_type=$4;
	context=$5;
	coverage=$6;
	printf "... Converting Cytosine Report to Bedgraph: %s\n" $output_type
	awk -v context=$context -v coverage=$coverage -f "$tool_dir""/"bismark2bedgraph.awk $report_infile > $bedgraph_outfile;
}
export -f awk_process;
# process function for the conversion from bedgraph to tdf
igv_process() {
	genome_chr_len=$1;
	bedgraph_infile=$2;
	tdf_outfile=$3;
	output_type=$4;
	printf "... Converting Bedgraph to TDF: %s\n" $output_type
	igvtools toTDF $bedgraph_infile $tdf_outfile $genome_chr_len > tmp_stdout;
}
export -f igv_process;

# parallelize the processes
parallel -k awk_process {1} {2} {3} {4} {5} {6} ::: "$tool_dir" :::+ "$infile_cov" ::: "${bedgraph_list[@]}" :::+ "${output_types[@]}" :::+ "${context_list[@]}" :::+ "${iscoverage[@]}";
parallel -k igv_process {1} {2} {3} {4} ::: "$igv_genome" ::: "${bedgraph_list[@]}" :::+ "${tdf_list[@]}" :::+ "${output_types[@]}";
