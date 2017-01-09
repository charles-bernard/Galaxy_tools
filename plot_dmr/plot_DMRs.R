#!/usr/bin/env Rscript

#############################################################################
# LOAD THE LIBRARIES
#############################################################################
library(optparse);
library(circlize);
library(data.table);
circos.clear()

#############################################################################
# READ THE OPTIONS 
#############################################################################
option_list = list(
	make_option("--context", dest="ctx", action="store"),
	make_option("--chr_len", dest="chr_len", action="store"),
	make_option("--hypo", dest="hypo", action="store"),
	make_option("--hyper", dest="hyper", action="store"),
	make_option("--output", dest="output", action="store")
);

#############################################################################
# EXTRACT THE OPTIONS 
#############################################################################
opt=parse_args(OptionParser(option_list=option_list));

#############################################################################
# MOHAMMED KASSIM'S SCRIPT
#############################################################################
if (opt$ctx=='CG'){
	hypercolor="red";
	hypocolor='red';
}
if (opt$ctx=='CHG'){
	hypercolor="darkblue";
	hypocolor='darkblue';
}
if (opt$ctx=='CHH'){
	hypercolor="darkgreen";
	hypocolor='darkgreen';
}
if (opt$ctx=='CDMR'){
	hypercolor="#202020";
	hypocolor='#202020';
}

datahyper <- fread(opt$hyper, header=TRUE, sep="\t");
datahypo <- fread(opt$hypo, header=TRUE, sep="\t");
lgChr <- fread(opt$chr_len, header=FALSE, sep="\t");

jpeg(opt$output, width=8000, height=4000, res=100, units="px", quality=100);

m <- matrix(c(1:5),5,1);
layout(m);
par(adj=0.5 ,bty="n", col.axis="dark blue", fg="dark blue", font.lab=50, ps=40, las=2, lab=c(20,10,10), mar=c(3,4,3,4));

for (c in 1:nrow(lgChr)){
	if ( lgChr[c,1] == "ChrC" || lgChr[c,1] == "ChrM" ) { 
		next 
	}
	lgchro=lgChr[c,2];
	datahyperchro=datahyper[datahyper$"#chr"==paste("Chr", c, sep=""),]
	datahypochro=datahypo[datahypo$"#chr"==paste("Chr", c, sep=""),]
	plot(c(1,lgchro), c(0,0), type="l", ylim=c(-2,2), ylab="", xlab="", xaxt="n", yaxt="n", xlim=c(0,30500000),
		main=paste("Chr",c,sep=""), lwd=4, cex.main=1.6)
	if (c==1){
		rect(0, -0.4, lgchro, 0.4, border="black", col="light gray")
		rect(11500000, -0.4, 17700000, 0.4, border="black", col="dimgray")
	}
	if (c==2){
		rect(0, -0.4, lgchro, 0.4, border="black", col="light gray")
		rect(1100000, -0.4, 7200000, 0.4, border="black", col="dimgray")
	}
	if (c==3){
		rect(0, -0.4, lgchro, 0.4, border="black", col="light gray")
		rect(10300000, -0.4, 17300000, 0.4, border="black", col="dimgray")
	}
	if (c==4){
		rect(0, -0.4, lgchro, 0.4, border="black", col="light gray")
		rect(1500000, -0.4, 2300000, 0.4, border="black", col="dimgray")
		rect(2800000, -0.4, 6300000, 0.4, border="black", col="dimgray")
	}
	if (c==5){
		rect(0, -0.4, lgchro, 0.4, border="black", col="light gray")
		rect(9000000, -0.4, 16000000, 0.4, border="black", col="dimgray")
	}
	if (length(datahyperchro[,1]) > 0 ){
		for (i in 1:length(datahyperchro[,1])){
			rect(datahyperchro[i,3], 0.8, datahyperchro[i,3], 1, border=hypercolor, lwd=8)
		}		 
	}
	if (length(datahypochro[,1]) > 0 ){
		for (i in 1:length(datahypochro[,1])){
			rect(datahypochro[i,3], -0.8, datahypochro[i,3], -1, border=hypocolor, lwd=8)
		}
	}
}
junk=dev.off()

# hypo_dmr <- fread(opt$hypo, header=TRUE, sep='\t');
# hyper_dmr <- fread(opt$hyper, header=TRUE, sep='\t');
# all_dmr <- rbind(hypo_dmr, hyper_dmr);
# colnames(all_dmr[,1:3])=c("chrom", "chromStart", "chromEnd");
# cytoband_table=data.table(all_dmr[,1:3], name=rep("", nrow(all_dmr)), gieStrain=rep("", nrow(all_dmr)));

# par(mar =c(1, 1, 1, 1));
# circos.initializeWithIdeogram(cytoband_table, sort.chr=TRUE);
# dmr_list=list(hypo_dmr[,1:3], hyper_dmr[,1:3]);
# circos.genomicRainfall(dmr_list, pch = 16, cex = 0.5, col =c( "red","green"));
# circos.genomicDensity(hypo_dmr[,1:3], col =c("red"), track.height = 0.08);
# circos.genomicDensity(hyper_dmr[,1:3], col =c("green"), track.height = 0.08);