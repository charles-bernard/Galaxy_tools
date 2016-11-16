#!/bin/bash

# read the options
TEMP=`getopt -o e::i:co: --long epi::,infile_cov:,context,output_dir:,tool_dir:,tdf,igv_genome: -n 'report' -- "$@"`
eval set -- "$TEMP"

# extract options and their arguments into variables.
while true ; do
	case "$1" in
		-e | --epi )
			#epi is the prefix of the output files names
			case "$2" in
				"" ) epi="current_job"; shift 2 ;;
				*) epi=$2; shift 2 ;;
			esac ;;
		-i | --infile_cov )
			#infile_cov is the filename of the cytosine report taken as input
			case "$2" in
				*) infile_cov=$2; shift 2 ;;
			esac ;;
		-c | --context ) 
			#context defines weither 2 or 4 bedgraphs are returned
			context=true; shift ;;
		--tdf )
			#tdf defines weither or not bedgraphs have to be converted into tdf files.
			tdf=true; shift ;;
		--igv_genome )
			#tdf conversion is achieved by igvtools and requires a file with the chrs_len of the genome
			case "$2" in
				*) igv_genome=$2; shift 2 ;;
			esac ;;
		-o | --output_dir )
			#output_dir in this galaxy tool is the tmp dir created by the wrapper.py
			case "$2" in
				*) output_dir=$2; shift 2 ;;
			esac ;;
		--tool_dir )
			#tool_dir is recquired to call other scripts stored in this directory
			case "$2" in
				*) tool_dir=$2; shift 2 ;;
			esac ;;
		-- ) shift; break ;;
		* ) echo "Internal error!"; exit 1 ;;
	esac
done

#IGV_path
#IGV_path="/users/biocomp/chbernar/galaxy_testing/database/dependencies/igvtools/2.3.32/geert-vandeweyer/package_igvtools_2_3_32/3c087cee3b8f/bin"

# define outputs according to options
if [[ "$context" = true ]]; then
	context_list=("CG" "CHG" "CHH")
	n="4"
	output_types=("CG" "CHG" "CHH" "coverage")
	bedgraph_list=("$output_dir""/""$epi""_CpG.bedgraph" "$output_dir""/""$epi""_CHG.bedgraph" "$output_dir""/""$epi""_CHH.bedgraph" "$output_dir""/""$epi""_coverage.bedgraph")
	if [[ "$tdf" = true ]]; then
		tdf_list=("$output_dir""/""$epi""_CpG.tdf" "$output_dir""/""$epi""_CHG.tdf" "$output_dir""/""$epi""_CHH.tdf" "$output_dir""/""$epi""_coverage.tdf")
	fi
else
	context_list=(".*")
	n="2"
	output_types=("CXX" "coverage")
	bedgraph_list=("$output_dir""/""$epi""_CXX.bedgraph" "$output_dir""/""$epi""_coverage.bedgraph")
	if [[ "$tdf" = true ]]; then
		tdf_list=("$output_dir""/""$epi""_CXX.tdf" "$output_dir""/""$epi""_coverage.tdf")
	fi
fi

# process
for (( i=0; i<$n; i++)); do
	printf "________________________________________________________________________\n"
	printf "Processing %s\n" ${output_types[$i]}
	printf "... Converting Cytosine Report to Bedgraph\n" 
	if (( i < n - 1 )); then #if not coverage:
		printf "#<Chr>\t<Start>\t<End>\t<Strand;Meth_ratio>\n" > "${bedgraph_list[$i]}"
		awk -v context="${context_list[$i]}" -v coverage="false" -f "$tool_dir"/bismark2bedgraph.awk $infile_cov >> "${bedgraph_list[$i]}"
	else
		printf "#<Chr>\t<Start>\t<End>\t<Coverage>\n" > "${bedgraph_list[$i]}"
		awk -v context="${context_list[$i]}" -v coverage="true" -f "$tool_dir"/bismark2bedgraph.awk $infile_cov >> "${bedgraph_list[$i]}"
	fi
	if [[ "$tdf" = true ]]; then
		printf "... Converting Bedgraph to Tdf\n"
		#"$IGV_path""/"igvtools toTDF "${bedgraph_list[$i]}" "${tdf_list[$i]}" "$igv_genome" > stdout_file  
		igvtools toTDF "${bedgraph_list[$i]}" "${tdf_list[$i]}" "$igv_genome" > stdout_file 
	fi
done