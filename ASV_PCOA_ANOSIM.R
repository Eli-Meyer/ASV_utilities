#######	Principle coordinate analaysis (PCoA) for ASV abundance data	#######
# Edit the input filename to match your filename, then run the script one section at a time.

# load libraries
library(vegan)
library(ape)

# read in your counts
filename<-"trimmed.ASVtable.txt"		# edit this to match your ASV table filename
counts<-read.table(filename, sep="\t", header=T)	
rownames(counts)<-counts$X
counts<-counts[,-1]
counts<-counts[,1:(ncol(counts)-5)]
head(counts[,1:4])				# view a subset of your input data
nrow(counts)					# number of ASVs in your data

# choose the number of ASVs to include in the analysis. e.g. 10 uses the
# top 10 most abundant ASVs. can be determined using "dada2 ASV barplot.R"
# or set manually.
asv.num<-19						# edit this to adjust ASV number
sym.counts<-counts[1:asv.num,]		
					
# log transform the counts data for pcoa
log.sym<-as.matrix(log10(sym.counts+1))	# using a log10 + 1 transformation
head(log.sym[,1:4])				# view a subset of your log-transformed data

# create a dissimilarity (distance) matrix
sym.dist<-vegdist(t(log.sym), method="manhattan")		# other methods available, see ?vegdist
tree<-hclust(sym.dist)
par(mar=c(12, 12, 5, 5))
# plot(as.dendrogram(tree), horiz=T, cex.wtf=0.2)		# optional - shows distance matrix as a cladogram

# establish a key
# edit as needed to reflect your sample codes
key<-as.vector(colnames(log.sym))
key<-gsub("FGB.+","TX",key,perl=TRUE)	# to assign all samples beginning with "FGB" to the 
key<-gsub("FK.+","FL",key,perl=TRUE)	# group "TX", use ["FGB.+","TX"]
key

# conduct statistical tests (analysis of similarity) to evaluate significance of
# differences in ASV profiles among groups
# p-value listed as "Significance" in the following output
sym.ano<-anosim(sym.dist, key)
summary(sym.ano)

# establish a color scheme for the PCoA plot. You may want to reflect significant differences 
# from anosim (above) in this color scheme. 
colvec<-key
# edit the sample labels and colors as appropriate for your data. e.g. to label all 
colvec<-gsub("TX", "red", colvec, perl=T)		# samples in group "FL" in green, set
colvec<-gsub("FL", "blue", colvec, perl=T)	# pattern to "FL" and replacement to "green"
t(cbind("sample"=colnames(counts), "color"=colvec))
# this is your color scheme

# run and plot the pcoa
sym.pcoa<-pcoa(sym.dist)$vectors
par(mar=c(12, 12, 5, 5))
par(las=1)
plot(sym.pcoa, ylab="Principal Coordinate 2", xlab="Principal Coordinate 1", 
	pch=19, col=colvec)
legend(legend=c("TX","FL"), col=c("red","blue"), pch=19, 
	horiz=T, xpd=T, x=-20, y=15)
