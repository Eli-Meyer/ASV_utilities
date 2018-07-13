####### A descriptive barplot showing ASV profiles by sample #######
# Edit the input filename to match your filename, then run the 
# script one section at a time.

# load libraries
# library(RColorBrewer)

# read in your counts
filename<-"trimmed.ASVtable.txt"		# edit this to match your ASV table filename
counts<-read.table(filename, sep="\t", header=T)	
rownames(counts)<-counts$X
counts<-counts[,-1]
head(counts[,1:4])	# view a subset of your input data
nrow(counts)		# number of ASVs in input

# convert symbiont sequence variant counts into proportions
sym.prop<-counts
for (a in 1:(ncol(counts)-5))
	{
	sym.prop[,a]<-sym.prop[,a]/sum(sym.prop[,a])
	}
head(sym.prop[,1:4])	# view a subset of your data as proportions

# choose the number of sequence variants to display. e.g. the minimum number 
# of ASVs that capture at least 95% of reads for all samples
minprop<-0.95		# choose your threshold here
for (a in 1:nrow(sym.prop))
	{
	cumsum<-(c(a, colSums(sym.prop[1:a,1:(ncol(sym.prop)-5)])))
	if (min(cumsum) >= minprop)
		{
		print(paste("Minimum number of ASVs that capture ",minprop*100,
			 "% of reads in all samples: ",a))
		break;
		}
	}
# the number of ASVs to be shown in the figure is automatically chosen based 
# on the analysis above. if you'd prefer to choose the number of ASVs manually, 
# you can adjust this here
anum<-a			# automatic default choice is "anum<-a". 
				# change it manually as e.g. "anum<-10"

# Establish a color scheme for the barplot. Here we encode sequence variants from a 
# "clade" (e.g. A, B, etc.) by different shades of a clade-specific color 
# (e.g. shades of blue for Bs). 
namevec<-sym.prop[1:anum,]$BestMatch
as.vector(namevec)		# the list of ASV IDs to be plotted
simpvec<-gsub("A.+","a",namevec,perl=T)
simpvec<-gsub("B.+","b",simpvec,perl=T)
simpvec<-gsub("C.+","c",simpvec,perl=T)
simpvec<-gsub("D.+","d",simpvec,perl=T)
colno<-max(table(simpvec))+1
yellows<-rainbow(n=colno, start=1/6, end=2/6)		
blues<-rainbow(n=colno, start=7/12, end=9/12)		
greens<-rainbow(n=colno, start=2/6, end=3/6, v=0.8)		
reds<-rainbow(n=colno, start=1/24, end=1/6)		
colvec<-simpvec
an<-bn<-cn<-dn<-0
for (c in 1:length(simpvec))
	{
	if (simpvec[c] == "a") {an<-an+1; colvec[c]<-yellows[an];}	# this section assigns clades
	if (simpvec[c] == "b") {bn<-bn+1; colvec[c]<-blues[bn];}	# to colors. Edit as needed
	if (simpvec[c] == "c") {cn<-cn+1; colvec[c]<-greens[cn];}	# for your dataset. e.g. switch
	if (simpvec[c] == "d") {dn<-dn+1; colvec[c]<-reds[dn];}	# "reds" and "greens".
	}


# make the bar plot
par(mar=c(12, 5, 5, 12))
par(las=1)
b<-barplot(as.matrix(sym.prop[1:anum,1:(ncol(sym.prop)-5)]), col=colvec,
	ylab="Proportion of sequence variants", names.arg=rep("", ncol(sym.prop)-5))
abline(h=0)
names<-colnames(sym.prop[,1:(ncol(sym.prop)-5)])
text(names, x=b, y=-0.1, xpd=T, srt=90, cex=0.8)

# make a key
legend(legend=namevec, fill=colvec, x=max(b)+1, y=1, xpd=T, horiz=F,
	cex=0.75, bty="n")
