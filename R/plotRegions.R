#!/usr/bin/env Rscript

####################################################################
# Visualization of sequence regions covered by RNA-RNA interactions
# predicted by IntaRNA.
#
# arguments: <IntaRNA-output-CSV> <1|2> <output-plot-file>
# 1 <IntaRNA-output-CSV> = ";"-separated CSV output of IntaRNA
# 2 <1|2> = suffix of "start,end,id" CSV cols to plot
# 3 <output-plot-file> = file name of the otuput figure suffixed
#     by one of ".pdf",".png",".svg",".eps",".ps",".jpeg",".tiff"
#
# example call:
#
# Rscript --slave -f plotRegions.R --args pred.csv 1 regions.png
#
# This script is part of the IntaRNA source code package. See
# respective licence and documentation for further information.
#
# https://github.com/BackofenLab/IntaRNA
#
####################################################################
# check and load dependencies
####################################################################

options(warn=-1)
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(ggalt))
suppressPackageStartupMessages(require(cowplot)) # cowplot starts with a note
options(warn=0)

theme_set(theme_cowplot())

####################################################################
# get command line arguments
####################################################################

args = commandArgs(trailingOnly=TRUE)
# check and parse
if (length(args)!=3) { stop("call with <intarna-csv-output> <1|2> <out-file-of-plot>", call.=FALSE) }

intarnaOutputFile = args[1];
if (!file.exists(intarnaOutputFile )) { stop("intarna-csv-output file '", intarnaOutputFile, "' does not exist!", call.=FALSE) }

seqNr = args[2];
if (seqNr != "1" &&  seqNr != "2") { stop("second call argument as to be '1' or '2' to specify which regions to plot"); }
# compile column names
id = paste("id",seqNr,sep="");
start = paste("start",seqNr,sep="");
end = paste("end",seqNr,sep="");

outFile = args[3];
fileExtensions = c(".pdf",".png",".svg",".eps",".ps",".jpeg",".tiff");
outFileExtOk = FALSE;
for( ext in fileExtensions ) { outFileExtOk = outFileExtOk || endsWith(outFile,ext); }
if ( !outFileExtOk ) {stop("<out-file-of-plot> has to have one of the following file extensions ",paste(fileExtensions,sep=" "), call.=FALSE);} 

# if set to some x-position, this will trigger the plotting of a vertical line at that location
xVline = NA;
#xVline = 250; # DEBUG TEST VALUE


####################################################################
# parse IntaRNA output
####################################################################

d = read.csv2( intarnaOutputFile )
# check if all columns present
for( x in c(id,start,end)) {
	if (!is.element(x, colnames(d))) { 
		stop("'",id,"' is not among the column names of '",intarnaOutputFile,"'", call.=FALSE); 
	}
}

####################################################################
# create count plot
####################################################################

allPos = c();
for( i in 1:nrow(d) ) {
	allPos = c( allPos, d[i,start]:d[i,end] );
}
allPos = as.data.frame(allPos,ncol=1)
#allPos # DEBUG OUT

coveragePlot =
		ggplot( allPos, aes(x=allPos, stat(count))) +
		geom_density() +
		ylab("coverage") +
		xlab("position") +
		scale_y_continuous(position = "right", expand=expand_scale(mult = c(0, .02))) +
		scale_x_continuous(expand = c(0, 0)) +
		theme(	axis.title.x=element_blank() )

if ( ! is.na(xVline)) {
	coveragePlot = coveragePlot +
			geom_vline(aes(xintercept=xVline));
}

####################################################################
# create region plot
####################################################################

dRegion = data.frame()
dRegion[1:nrow(d),1:3] = d[,c(start,end,id)]
dRegion[1:nrow(d),4] = factor(sprintf("%08d",nrow(d):1))
colnames(dRegion) = c("start","end","id","idx");
#dRegion # DEBUG OUT

yLabelScale = 0.6 # if you have to alter for more/less sequence IDs per inch; see <plotHeight> below

regionPlot = 
		ggplot(dRegion, aes(x=start,xend=end,y=idx)) +
		geom_dumbbell(color="dodgerblue", size=2) +
		xlab("position") + 
		scale_x_continuous( expand = c(0, 0) ) +
		ylab("") +
		scale_y_discrete(position = "right", breaks=dRegion$idx, labels=dRegion$id) +
		geom_vline(aes(xintercept=min(allPos))) +
		theme(panel.grid.major.y=element_line(size=0.7,color="lightgray")
				, axis.text.y=element_text(size=rel(yLabelScale))
		)

if ( ! is.na(xVline)) {
	regionPlot = regionPlot + 
			geom_vline(aes(xintercept=xVline));
}

####################################################################
# plot to file
####################################################################

plotWidth = 6
plotHeightDensity = 2
plotHeight = plotHeightDensity + max(2,nrow(d)/9)
plotHeightDensityRel = plotHeightDensity / plotHeight


plot_grid( coveragePlot, regionPlot
		, nrow=2, ncol=1
		, align = "hv"
		, axis = "r"
		, rel_heights= c( plotHeightDensityRel, 1.0-plotHeightDensityRel)
)

ggsave( outFile
		, width= plotWidth
		, height= plotHeight
);

#############################################################EOF
