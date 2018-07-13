#!/usr/bin/env perl
# written by E Meyer, eli.meyer@science.oregonstate.edu
# distributed without any guarantees or restrictions

# -- check arguments and print usage statement
$scriptname=$0; $scriptname =~ s/.+\///g;
$usage = <<USAGE;
Draws the specified number of sequences randomly from a pair of input files.
This provides a randomly selected set of paired reads. 
Usage: $scriptname -f f.input -r r.input -n num_seq -o f.output -p r.output
Required arguments:
	-f f.input	name of the "R1" input file (forward reads) from which sequences will be randomly drawn.
	-r r.input	name of the corresponding "R2" input file (reverse reads).
	-n num_seq	number of sequences to draw
	-o f.output	a name for the R1 (forward) output file
	-p r.output	a name for the R1 (forward) output file
USAGE
if ($#ARGV < 4 || $ARGV[0] eq "-h") {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# -- module and executable dependencies
$mod1="Getopt::Std";
unless(eval("require $mod1")) {print "$mod1 not found. Exiting\n"; exit;}
use Getopt::Std;

# get variables from input
getopts('f:r:n:o:p:h');	# in this example a is required, b is optional, h is help
if (!$opt_f || !$opt_r || !$opt_n || !$opt_o || !$opt_p || $opt_h) {print "\n", "-"x60, "\n", $scriptname, "\n", $usage, "-"x60, "\n\n"; exit;}

# name variables from user input
$infile = $opt_f;
$revfile = $opt_r;
$todraw = $opt_n;
$outfile1 = $opt_o;
$outfile2 = $opt_p;

# count number of sequences in input
$nolines = `wc -l $infile`;
$nolines =~ s/\s.+//g;
$revlines = `wc -l $revfile`;
$revlines =~ s/\s.+//g;
if ($nolines ne $revlines)
	{
	print "$infile and $revfile have different numbers of sequences. Check input and try again.\n";
	exit;
	}
chomp($nolines);
$noseqs = $nolines / 4;

if ($todraw >= $noseqs)
	{
	print "\nThere are only $noseqs sequences in $infile and you've asked to randomly draw $todraw.\n\n";
	exit;
	}

# generate list of sequence numbers to extract
%rh = ();
for ($a=1; $a<=$todraw; $a++)
	{
	$randi = 0;
	while ($randi eq 0 || exists($rh{$randi}))
		{
		$randi = int(rand($noseqs+1));
		}
	$rh{$randi}++;
#	print $randi, "\n";
	}

# loop through R1 fastq file and print out randomly chosen sequence numbers
open (IN, $infile);
open (OUT, ">$outfile1");
$state = 0; $seqcount = 0; $count = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) 
		{
		$seqcount++;
		if (exists($rh{$seqcount}))
			{$state = 1;}
		else 	{$state = 0;}
		}
	if ($state eq 1) {print OUT $_, "\n";}
	}
close(IN);

# loop through R1 fastq file and print out randomly chosen sequence numbers
open (IN, $revfile);
open (OUT, ">$outfile2");
$state = 0; $seqcount = 0; $count = 0;
while(<IN>)
	{
	chomp;
	$count++;
	if ($count==1) {$ss = substr($_, 0, 4);}
	if ($_ =~ /^$ss/) 
		{
		$seqcount++;
		if (exists($rh{$seqcount}))
			{$state = 1;}
		else 	{$state = 0;}
		}
	if ($state eq 1) {print OUT $_, "\n";}
	}
close(IN);
