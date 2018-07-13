A collection of scripts for ITS amplicon analysis.

-------------------------
blastnID.pl
-------------------------

------------------------------------------------------------
blastnID.pl
Compares a set of DNA sequences with a curated DNA database to infer the biological origin of
each sequence. 
Usage: blastnID.pl -q queries -d DB -s min_score -i input.txt -o output.txt
Required arguments:
	-q queries	the file of sequences to be identified, FASTA format
	-d DB		a BLAST formatted nucleotide database
	-s min_score	alignment threshold (bit-score), higher scores reflect greater similarity. 
	-i input.txt	name of the input ASV table (from dada2.R).
	-o output.txt	a name for the annotated ASV table output (tab-delimited text file).
------------------------------------------------------------

-------------------------
FindPairedReads.pl
-------------------------

------------------------------------------------------------
FindPairedReads.pl
Finds valid paired reads in a pair of files supplied by the user (FASTQ). 
Paired reads are written in the same order in a pair of files named
after the input files, with a prefix of 'paired.'
If the option "b" is chosen, orphans are also written to files named
after the input, prefixed with 'orphan.'
Usage: FindPairedReads.pl -f F_reads -r R_reads <OPTIONS>
Required arguments:
	-f F_reads	input, forward reads (R1)(FASTQ)
	-r R_reads	input, reverse reads (R2)(FASTQ)
Options:
	-o option	choose 'p' to print only paired reads (default),
			or 'b' to print both paired and orphan reads.
	-g F_pair	a name for the paired forward (R1) reads output file (default: paired.input)	
	-i R_pair	a name for the paired reverse (R2) reads output file (default: paired.input)	
	-j F_orph	a name for the orphan forward (R1) reads output file (default: orphan.input)	
	-k R_orph	a name for the orphan reverse (R2) reads output file (default: orphan.input)	
------------------------------------------------------------

-------------------------
QualFilterFastq.pl
-------------------------

------------------------------------------------------------
QualFilterFastq.pl

Removes reads containing too many low quality basecalls from a set of short sequences
Output:  high-quality reads in FASTQ format
Usage:   QualFilterFastq.pl -i input -m min_score -x max_LQ -o output
Required arguments:
	-i input	raw input reads in FASTQ format
	-m min_score	quality scores below this are considered low quality (LQ)
	-x max_LQ	reads with more than this many LQ bases are excluded
	-o output	name for ourput file of HQ reads in FASTQ format
------------------------------------------------------------

-------------------------
QualViewFastq.pl
-------------------------

------------------------------------------------------------
QualViewFastq.pl
Shows descriptive statistics for quality scores by position in a file of 
DNA sequences provided by the user. The output is a box plot in text format,
following the same conventions as the html output from FASTQC. 
Usage:   QualViewFastq.pl -i reads <OPTIONS> 
Required arguments:
	-i reads	the FASTQ file to be analyzed (the script expects reads of equal length)

Options:
	-w <NUMBER>	window	size of the region shown (basepairs); longer reads will be shown in 
			blocks of this size (default = 100)
------------------------------------------------------------

-------------------------
random_paired_fastq.pl
-------------------------

------------------------------------------------------------
random_paired_fastq.pl
Draws the specified number of sequences randomly from a pair of input files.
This provides a randomly selected set of paired reads. 
Usage: random_paired_fastq.pl -f f.input -r r.input -n num_seq -o f.output -p r.output
Required arguments:
	-f f.input	name of the "R1" input file (forward reads) from which sequences will be randomly drawn.
	-r r.input	name of the corresponding "R2" input file (reverse reads).
	-n num_seq	number of sequences to draw
	-o f.output	a name for the R1 (forward) output file
	-p r.output	a name for the R1 (forward) output file
------------------------------------------------------------

-------------------------
ScreenFastqByPrimer.pl
-------------------------

------------------------------------------------------------
ScreenFastqByPrimer.pl
Filters a set of short reads in FASTQ format, excluding any reads that lack the 
expected sequence at the 5' end of the read. This is typically used to identify
valid amplicons in amplicon sequence analysis. 
Usage: ScreenFastqByPrimer.pl -i input -s sequence -m mismatches <OPTIONS>
Required arguments:
	-i input	name of the input file, FASTQ
	-s sequence	expected sequence (primer sequence)
	-m mismatches	maximum number of mismatches allowed between expected and observed sequences	
Options:
	-p pass.name	a name for the output file of sequences passing this filter (default: target.input)
	-f fail.name	a name for the output file of sequences failing this filter (default: nontarget.input)
------------------------------------------------------------

-------------------------
TruncateFastq.pl
-------------------------

------------------------------------------------------------
TruncateFastq.pl
Truncates a set of short reads in FASTQ format to keep the region specified
Usage:   TruncateFastq.pl -i input -s start -e end -o output
Required arguments:
	-i input	file of short reads to be filtered, fastq format
	-s start	beginning of the region to keep, nucleotide position
	-e end		end of the region to keep, nucleotide position
	-o output	a name for the output file (fastq format)
------------------------------------------------------------

