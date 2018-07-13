#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Compares a set of DNA sequences with a curated DNA database to infer the biological origin of
each sequence. 
Usage: $scriptname -q queries -d DB -s min_score -i input.txt -o output.txt
Required arguments:
	-q queries	the file of sequences to be identified, FASTA format
	-d DB		a BLAST formatted nucleotide database
	-s min_score	alignment threshold (bit-score), higher scores reflect greater similarity. 
	-i input.txt	name of the input ASV table (from dada2.R).
	-o output.txt	a name for the annotated ASV table output (tab-delimited text file).
USAGE
if ($#ARGV < 4 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="File::Which";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use File::Which;
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;
$mod1="Bio::SearchIO";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Bio::SearchIO;

# use this block and edit to check for executables required by the script
$dep1 = "blastn";
unless (defined(which($dep1))) {print $dep1, " not found. Exiting.\n"; exit;}

# get variables from input
getopts('q:d:s:i:o::h');	# in this example a is required, b is optional, h is help
if (!$opt_q || !$opt_d || !$opt_s || !$opt_i || !$opt_o || $opt_h) 
	{print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}
$qfile = $opt_q;
$dbfile = $opt_d;
$critscore = $opt_s;
$tabfile = $opt_i;
$outfile = $opt_o;

# usage
use Bio::SearchIO;

# run blast
if (-e "$qfile.blastn.txt")
	{
	print "$qfile.blastn.txt already exists. Using that file...\n";
	}
else
	{
	print "Beginning BLASTN search...";
	system("date");
	system("blastn -db $dbfile -query $qfile -num_alignments 10 -num_descriptions 10 -out $qfile.blastn.txt");
	print "Finished BLASTN search.\n";
	}

# parse blast report
print "Parsing BLASTN report...\n";
system("date");
$report = new Bio::SearchIO(-file=>"$qfile.blastn.txt", -format=>"blast");
while ($result = $report->next_result)
	{
	$qid = $result->query_name;
	$qid =~ s/\;.+//;
	$hid = $score = $eval = $besthit = $bestscore = $bestE = $nexthit = $nextscore = $nextE = $erat = "";
	$hitcount = 0;
	$queries++;
	while ($hit = $result->next_hit)
		{
		$hid = $hit->accession;
		$eval = $hit->significance;
		$score = $hit->bits;
		if ($score >= $critscore)
			{
			if ($hitcount == 0)
				{
				$besthit = $hid;
				$bestscore = $score;
				$bestE = $eval;
				}
			elsif ($hitcount == 1)
				{
				$nexthit = $hid;
				$nextscore = $score;
				$nextE = $eval;
				}
			$hitcount++;
			}
		}
	if ($nexthit ne "")
		{
		$erat = $bestE/$nextE;
		}
	else	{$erat = 0;}
#	print join ("\t", ($qid, $besthit, $bestscore, $erat, $nexthit, $nextscore)), "\n";
	if ($besthit ne "") 
		{
		$sigmatch++;
		$bh{$qid}{"besthit"} = $besthit;
		$bh{$qid}{"bestscore"} = $bestscore;
		$bh{$qid}{"eratio"} = $erat;
		$bh{$qid}{"nexthit"} = $nexthit;
		$bh{$qid}{"nextscore"} = $nextscore;
		}
	}
print "Finished parsing BLASTN report.\n";
system("date");

# read OTU table
system("dos2unix $tabfile");
open(TAB,$tabfile);
open(OUT, ">$outfile");
while(<TAB>)
	{
	chomp;
	$rownum++;
# write annotated OTU table
	if ($rownum == 1) 
		{
		print OUT $_, "\t";
		print OUT join("\t", qw{BestMatch BestScore EvalRatio NextHit NextScore}), "\n";
		next;
		}
	@cols = split("\t", $_);
	$id = $cols[0];
	print OUT $_, "\t";
	if ($bh{$id}{"besthit"} ne "")
		{
		print OUT join ("\t", ($bh{$id}{"besthit"}, $bh{$id}{"bestscore"}, $bh{$id}{"eratio"}, 
			$bh{$id}{"nexthit"}, $bh{$id}{"nextscore"})), "\n";
		}
	else
		{
		print OUT join ("\t", qw{- - - - -}), "\n";
		}
	}	


# write summary 
print "$queries sequences found in $qfile\n";
print "$sigmatch of these matched one or more record in $dbfile with scores >= $critscore\n";
