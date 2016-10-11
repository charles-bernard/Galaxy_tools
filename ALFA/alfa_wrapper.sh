#!usr/bin/bash 

#########################################################################################################
# ARGUMENTS FROM alfa_wrapper.xml                                                                       #
#########################################################################################################
galaxyRoot=$1;
toolDir=$2
configFile=$3;
logReport=$4;
sed -i -e '/^$/d; s/\t//g;' $configFile;
printf "__________________________________________________________________\n\n" > $logReport
printf "                          ALFA CONFIG                             \n" >> $logReport
printf "__________________________________________________________________\n" >> $logReport
cat $configFile >> $logReport

#########################################################################################################
# INITIALIZATION OF THE VARIABLES from $configFile                                                      #
#########################################################################################################
#_INPUT1
annotationSource=`grep -P '^annotationSource ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
if [ "$annotationSource" == "personal_gtf" ]; then
	annotationFile=`grep -P '^annotationFile ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
elif [ "$annotationSource" == "built_in_index" ]; then
	built_in_index_prefix=`grep -P '^built_in_index_prefix ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
else
	strandedIndex=`grep -P '^strandedIndex ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
	unstrandedIndex=`grep -P '^unstrandedIndex ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
fi

#_INPUT2
readsType=`grep -P '^readsType ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
readsFileList=`grep -P '^readsFile\[[0-9]+\] ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
readsLabelList=`grep -P '^readsLabel\[[0-9]+\] ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;

#_OUTPUT CHOICES
plotChoice=`grep -P '^plotChoice ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
countFileChoice=`grep -P '^countFileChoice ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
indexChoice=`grep -P '^indexChoice ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;

#_OUTPUT OPTIONS
strandness=`grep -P '^strandness ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
categoriesDepth=`grep -P '^categoriesDepth ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
plotFormat=`grep -P '^plotFormat ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
plotThresholdChoice=`grep -P '^plotThresholdChoice ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
if [ "$plotThresholdChoice" == "True" ]; then
	yMin=`grep -P '^yMin ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
 	yMax=`grep -P '^yMax ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
fi

#_OUTPUT FILES
if [ "$plotChoice" == "True" ]; then 
	if [ "$plotFormat" == "pdf" ]; then
		outputPdf=`grep -P '^outputPdf ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
	elif [ "$plotFormat" == "svg" ]; then
		outputCategoriesSvg=`grep -P '^outputCategoriesSvg ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
		outputBiotypesSvg=`grep -P '^outputBiotypesSvg ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
	else
		outputCategoriesPng=`grep -P '^outputCategoriesPng ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
		outputBiotypesPng=`grep -P '^outputBiotypesPng ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
	fi
fi
if [ "$countFileChoice" == "True" ]; then 
	outputCountFile=`grep -P '^outputCountFile ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
fi
if [ "$indexChoice" == "True" ]; then 
	outputStrandedIndex=`grep -P '^outputStrandedIndex ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
	outputUnstrandedIndex=`grep -P '^outputUnstrandedIndex ?=' $configFile | awk 'BEGIN{FS="= ?"} {print $2}'`;
fi

#########################################################################################################
# CREATION OF A TMP DIRECTORY FOR THE OUTPUT FILES OF ALFA AND cd                                       #
#########################################################################################################
outputDirectory=`mktemp -d "$galaxyRoot"/database/tmp/tmpXXXXXX`;
if [ -d $outputDirectory ]; then
	chmod -R ugo+wrx $outputDirectory;
	rm -R $outputDirectory;
fi
mkdir $outputDirectory;
chmod -R ugo+wrx $outputDirectory;
cd $outputDirectory;

#########################################################################################################
# TEST OF INPUT1                                                                                        #
#########################################################################################################
if [ "$annotationSource" == "index" ]; then
	#need to copy the files.dat to .*index because ALFA requires the extension ".(un)stranded.index"
	index="index"
	cp $strandedIndex $index".stranded.index"
	cp $unstrandedIndex $index".unstranded.index"
fi

#########################################################################################################
# TEST OF INPUT2 AND DETERMINATION OF PYTHON READS INPUT ARGUMENT                                       #
#########################################################################################################
readsListLen=`echo "$readsFileList" | wc -l`;
readsInput="";
for (( i = 1; i <= readsListLen; i++ )) do
	readsFile[$i]=`echo "$readsFileList" | awk -v i=$i 'NR==i'`;
	readsLabel[$i]=`echo "$readsLabelList" | awk -v i=$i 'NR==i' | sed -e 's/ /_/g'`;
	if [ "$readsType" == "bam" ]; then
		bamSorted=`samtools view -H "${readsFile[$i]}" | grep -c 'SO:unsorted'`
		if [ "$bamSorted" != "0" ] ; then
			samtools sort ${readsFile[$i]} ${readsFile[$i]}
		fi
	else
		#need to copy the file.dat to tmp.bedgraph because ALFA requires the extension ".bedgraph"
		bedgraphFile="tmpBedgraph_"$i
		cp ${readsFile[$i]} $bedgraphFile".bedgraph"
		readsFile[$i]=$bedgraphFile
	fi
	if [ "${readsLabel[$i]}" == "" ]; then
		readsLabel[$i]="sample_""$i";
	fi
	readsInput=$readsInput" "${readsFile[$i]}" "${readsLabel[$i]};
done

#########################################################################################################
# DETERMINATION OF THE APPROPRIATE SCRIPTS ARGUMENTS                                                    #
#########################################################################################################
scriptPath="$toolDir/";
if [ "$annotationSource" == "index" ]; then
	scriptInput="-g $index -i ""$readsInput";
elif [ "$annotationSource" == "built_in_index" ]; then
	scriptInput="-g $built_in_index_prefix -i ""$readsInput";
else
	scriptInput="-a $annotationFile -i ""$readsInput";
fi
if [ "$readsType" = "bedgraph" ]; then
	scriptInput=$scriptInput" --bedgraph";
fi
scriptStrandness="-s "$strandness
scriptCategoriesDepth="-d "$categoriesDepth
if [ "$plotChoice" == "True" ]; then
	if [ "$plotFormat" == "pdf" ]; then
		scriptPlotOutput="--pdf plotFile.pdf";
	else
		scriptPlotOutput="--"$plotFormat" plotFile";
	fi
	if [ "$plotThresholdChoice" == "True" ]; then
		scriptPlotOutput=$scriptPlotOutput" -t ""$yMin"" ""$yMax"
	fi
else
	scriptPlotOutput="--n";
fi

#########################################################################################################
# DISPLAY ALFA PROCESS                                                                                  #
#########################################################################################################
printf "__________________________________________________________________\n\n" >> $logReport
printf "                          ALFA PROCESS                            \n" >> $logReport
printf "__________________________________________________________________\n" >> $logReport

if [ "$plotChoice" == "False" ] && [ "$countFileChoice" == "False" ] && [ "$indexChoice" == "False" ]; then
cat <<error 1>&2

No output to return. 
Process Aborted
error
exit 0
fi

printf "Command:\n" >> $logReport
echo "python ""$scriptPath"ALFA.py $scriptInput $scriptStrandness $scriptCategoriesDepth $scriptPlotOutput >> $logReport;
printf "\n******************************************************************\n" >> $logReport
printf "Temporary Output Directory:\n" >> $logReport
echo $outputDirectory >> $logReport
printf "\n******************************************************************\n" >> $logReport
printf "ALFA prompt:\n" >> $logReport
python "$scriptPath"ALFA.py $scriptInput $scriptStrandness $scriptCategoriesDepth $scriptPlotOutput >> $logReport 2>errorFile;
printf "\n******************************************************************\n" >> $logReport

#########################################################################################################
# REDIRECTION OF ERRORS - TMP SOURCE ALFA.PY MUST BE CORRECTED SOON                                     #
#########################################################################################################
if [[ -s errorFile ]]; then
	#When the option --n is enabled, alfa prints '### End of the program' in stderr even if the process worked-
	#The following lines to avoid the tool from crashing in this case
	endProgram=`grep -c '### End of program' errorFile`
	if [ "$endProgram" == "0" ]; then
		#When alfa prints '### End of program' in stdout, all the messages in stderr are considered
		#as warnings and not as errors. True errors make the script exits with code "2"
		endProgram=`grep -c '### End of program' $logReport`
		if [ "$endProgram" == "0" ]; then
 			>&2 printf "The script ALFA.py encountered the following error:\n\n"
			>&2 cat errorFile
			printf "ALFA error:\n" >> $logReport
			cat errorFile >> $logReport
			printf "\n******************************************************************\n" >> $logReport
 			exit 2
 		else
 			>&2 printf "The script ALFA.py encountered the following warning:\n\n"
 			>&2 cat errorFile 
 			printf "ALFA warning:\n" >> $logReport
 			cat errorFile >> $logReport
			printf "\n******************************************************************\n" >> $logReport
 		fi
 	fi
fi

#########################################################################################################
# OUTPUT REDIRECTIONS                                                                                   #
#########################################################################################################
if [ "$plotChoice" == "True" ]; then
	if [ "$plotFormat" == "pdf" ]; then
		mv "plotFile.pdf" $outputPdf;
	elif [ "$plotFormat" == "png" ]; then
		mv "plotFile.categories.png" $outputCategoriesPng;
		mv "plotFile.biotypes.png" $outputBiotypesPng;
	else 
		mv "plotFile.categories.svg" $outputCategoriesSvg;
		mv "plotFile.biotypes.svg" $outputBiotypesSvg;
	fi
fi
if [ "$countFileChoice" == "True" ]; then
	> countFile;
	for (( i = 1; i <= readsListLen; i++ )) do
		printf "##LABEL: "${readsLabel[$i]}"\n\n" >> countFile;
		cat ${readsLabel[$i]}".categories_counts" >> countFile;
		printf "__________________________________________________________________\n" >> countFile;
	done
	mv countFile $outputCountFile;
fi
if [ "$indexChoice" == "True" ]; then
	if [ "$annotationSource" == "index" ]; then
		mv $strandedIndex $outputStrandedIndex
		mv $unstrandedIndex $outputUnstrandedIndex
	elif [ "$annotationSource" == "built_in_index" ]; then
		cp $built_in_index_prefix".stranded.index" $outputStrandedIndex
		cp $built_in_index_prefix".unstranded.index" $outputUnstrandedIndex
	else
		annotationFileName=`grep -P -o '[^/]*\.dat$' <<< $annotationFile`
		mv $annotationFileName".stranded.index" $outputStrandedIndex
		mv $annotationFileName".unstranded.index" $outputUnstrandedIndex
	fi
fi