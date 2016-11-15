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

#IGV_path
#IGV_path="/users/biocomp/chbernar/galaxy_testing/database/dependencies/igvtools/2.3.32/geert-vandeweyer/package_igvtools_2_3_32/3c087cee3b8f/bin"

# define outputs according to options
if [[ "$context" = true ]]; then
	context_list=("CG" "CHG" "CHH")
	output_types=("CG" "CHG" "CHH" "coverage")
	bedgraph_list=("$output_dir""/""$epi""_CpG.bedgraph" "$output_dir""/""$epi""_CHG.bedgraph" "$output_dir""/""$epi""_CHH.bedgraph" "$output_dir""/""$epi""_coverage.bedgraph")
	tdf_list=("$output_dir""/""$epi""_CpG.tdf" "$output_dir""/""$epi""_CHG.tdf" "$output_dir""/""$epi""_CHH.tdf" "$output_dir""/""$epi""_coverage.tdf")
	n="4"
else
	context_list=(".*")
	output_types=("CXX" "coverage")
	bedgraph_list=("$output_dir""/""$epi""_CXX.bedgraph" "$output_dir""/""$epi""_coverage.bedgraph")
	tdf_list=("$output_dir""/""$epi""_CXX.tdf" "$output_dir""/""$epi""_coverage.tdf")
	n="2"
fi

# process
for (( i=0; i<$n; i++)); do
	printf "________________________________________________________________________\n"
	printf "Processing %s\n" ${output_types[$i]}
	printf "... Converting Cytosine Report to Bedgraph\n" 
	if (( i < n - 1 )); then
		#if not coverage:
		#printf "track type=bedGraph name=%s Coverage description=%s Coverage\n" "$epi""_""${context_list[$i]}" "$epi""_""${context_list[$i]}" > "${bedgraph_list[$i]}"
		printf "#<Chr>\t<Start>\t<End>\t<Strand;Meth_ratio>\n" > "${bedgraph_list[$i]}"
		awk -v context="${context_list[$i]}" -v coverage="false" -f "$tool_dir"/bismark2bedgraph.awk $infile_cov >> "${bedgraph_list[$i]}"
	else
		#printf "track type=bedGraph name=%s Coverage description=%s Coverage\n" "$epi""_""${context_list[$i]}" "$epi""_""${context_list[$i]}" > "${bedgraph_list[$i]}"
		printf "#<Chr>\t<Start>\t<End>\t<Coverage>\n" > "${bedgraph_list[$i]}"
		awk -v context="${context_list[$i]}" -v coverage="true" -f "$tool_dir"/bismark2bedgraph.awk $infile_cov >> "${bedgraph_list[$i]}"
	fi
	if [[ "$tdf" = true ]]; then
		printf "... Converting Bedgraph to Tdf\n"
		#"$IGV_path""/"igvtools toTDF "${bedgraph_list[$i]}" "${tdf_list[$i]}" "$igv_genome" > stdout_file  
		igvtools toTDF "${bedgraph_list[$i]}" "${tdf_list[$i]}" "$igv_genome" > stdout_file 
	fi
done
