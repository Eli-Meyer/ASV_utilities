####### Testing for differences in abundance of specific ASVs #######
# Edit the input filename to match your filename, then run the script one section at a time.

# choose your significance threshold (FDR)
p.thd<-0.1

# load libraries
library(DESeq2)
library(gplots)

# read in your counts
filename<-"trimmed.ASVtable.txt"		# edit this to match your ASV table filename
counts<-read.table(filename, sep="\t", header=T)	
rownames(counts)<-paste(counts$X," (",counts$BestMatch,")",sep="")
counts<-counts[,2:(ncol(counts)-5)]
head(counts[,1:4])				# view a subset of your input data
nrow(counts)					# number of ASVs in your data

# choose the number of ASVs to include in the analysis. e.g. 10 uses the
# top 10 most abundant ASVs. can be determined using "dada2 ASV barplot.R"
# or set manually.
asv.num<-19						# edit this to adjust ASV number
sym.counts<-counts[1:asv.num,]	
nrow(sym.counts)					# number of ASVs remaining

# (optional) Exclude ASVs present at very low abundance from your analysis
# These are unlikely to reach significance but incur penalties in MTC
# here we exclude any ASVs present at < 10 reads per sample on average
covthd<-2						# edit this to set minimum coverage threshold
sym.counts<-sym.counts[rowSums(sym.counts)/ncol(sym.counts)>=covthd,]
nrow(sym.counts)					# number of ASVs remaining

# make a key describing your samples
samples<-as.vector(colnames(counts))
origins<-gsub("FGB.+","TX",samples,perl=TRUE)	# to assign all samples beginning with "FGB" to the 
origins<-gsub("FK.+","FL",origins,perl=TRUE)	# group "TX", use ["FGB.+","TX"]
key<-data.frame(cbind("sample"=samples, "origin"=origins))
head(key)

# construct a DESeq expression object
# this is where you define the complete model
deobj<-DESeqDataSetFromMatrix(sym.counts, key, design =~ origin)

# estimate size factors
deobj<-estimateSizeFactors(deobj, type="poscounts")

# estimate dispersions
deobj<-estimateDispersions(deobj)

vsd <- assay(varianceStabilizingTransformation(deobj, fitType="local"))
head(vsd)
nrow(vsd)

# compare complete model with reduced model to test significance of the term removed
fit1<-nbinomLRT(deobj, reduced =~ 1)
res1<-results(fit1, alpha=p.thd, cooksCutoff=F)		# turned off outlier filtering
summary(res1)
# this is a summary of the statistical tests
nres1<-res1[!is.na(res1$padj),]
sig1<-nres1[nres1$padj<p.thd,]
# view a subset of the ASVs present at significantly different proportions between groups
# number of significant ASVs:
nrow(sig1)
head(sig1)

# heatmap of ASVs present at significantly difference proportions between groups
# this color scheme reflects raw abundance
par(mar=c(5,5,12,12))
sig.vsd<-vsd[rownames(vsd) %in% rownames(sig1),]
rownames(sig.vsd)
library("RColorBrewer")
hmcol<-colorRampPalette(c("black","black", "black", "red3", "yellow")) (100)
heatmap.2(sig.vsd, col=hmcol, scale="none", margins=c(12,12), cexRow=1,
	trace="none", key.title=NA, keysize=1, density.info="none",
	key.xlab="log10(abundance+1)", key.par=list(mar=c(4,2,4,0)))
