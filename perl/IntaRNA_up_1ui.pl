#!/usr/bin/env/perl -w

use strict;
use Getopt::Std;
use File::Basename qw( dirname );

my $intaRNAbinPath = dirname(__FILE__)."/";

my %args;

getopt("t:m:v:s:p:w:L:l:T:u:", \%args);

if (defined $args{h} && $args{h}==1) {
print "
\n
Synopsis    : IntaRNA [-t fasta_file] [-m fasta_file] "#[-v energy] 
."[-o]\n
                      [-s number] [-n]  [-h] "#[-u[1|2] number] 
."[-u number]\n
                      [-p number] "#[-f pos,pos] 
."[-T temp] "#[-U] [-P]
."[-w size]\n
                      [-L distance] [-l length] \n"
#[-a weight] [-b weight]\n
#                      [-c threshold] target-RNA-seq binding-RNA-seq\n
."\n
Description : User interface wrapper for intaRNA v1 like calls.\n
\n
Options     :\n
 === General parameters ===\n
              -t fasta_file    : use fasta file of target sequences\n
              -m fasta_file    : use fasta file of binding sequences\n
              -v energy        : outputs all results below energy in kcal/mol\n
              -o               : detailed output\n
              -s number        : max. number of calculated suboptimal results\n
                                 (default:0)\n
              -n               : use no seed constraint\n
              -h               : this help\n
\n
 === Seed parameters ===\n
              -p number        : exact number of paired bases in the seed region\n
                                 (default:6)\n"
#              -u[1|2] number   : max. number of unpaired bases in the seed\n
#                                 region in\n
#                                 1: the first  sequence (default:0)\n
#                                 2: the second sequence (default:0)\n
."               -u number        : max. number of unpaired bases in the seed\n
                                 region in both sequences (default:0)\n"
#              -f number,number : search for seed in binding RNA (e.g. ncRNA)\n
#                                 in region between positions start,end\n
#                                 (given in 5' to 3' direction counting from 1)\n

."\n
 === RNA folding parameters ===\n
              -T temp          : temperature in Celsius (default: 37°C)\n"
#              -U               : use RNAup to compute ED values of binding RNA\n
#                                 (default)\n
#              -P               : use RNAplfold to compute ED values of target RNA\n
#                                 (default)\n
."              -w size          : window size for computation of ED values\n
                                 with RNAplfold (default: length of target RNA)\n
              -L distance      : max. distance of two paired bases for\n
                                 computation of ED values with RNAplfold\n
                                 (default: window size)\n
              -l length        : max. length of hybridized region, mainly used\n
                                 for efficient computation (default: window size)\n"
#              -a weight        : weight for ED values of target RNA in energy\n
#                                 (default: 1.0)\n
#              -b weight        : weight for ED values of binding RNA in energy\n
#                                 (default: 1.0)\n
#              -c threshold     : threshold for seed accessibility, requires u=0\n
#                                 EXPERIMENTAL FEATURE, (default: -1.0)\n
."";
exit 0;
}


# generate intaRNA 2 call
my $intaRNA2call = "";
if (defined $args{t}) { $intaRNA2call .= " -t ".($args{t}); }
if (defined $args{m}) { $intaRNA2call .= " -q ".($args{m}); }
if (defined $args{v}) { $intaRNA2call .= " --outMaxE=".($args{v}); }
if (defined $args{s}) { $intaRNA2call .= " -n ".($args{s}+1); }
if (defined $args{p}) { $intaRNA2call .= " --seedBP=".($args{p}); }
if (defined $args{T}) { $intaRNA2call .= " --temperature=".($args{T}); }
if (defined $args{w}) { 
	$intaRNA2call .= " --tAccW=".($args{w}); 
} else {
	$intaRNA2call .= " --tAccW=0"; 
}
# always full sequence length for query
$intaRNA2call .= " --qAccW=0"; 
if (defined $args{L}) { 
	$intaRNA2call .= " --tAccL=".($args{L}); 
} else {
	$intaRNA2call .= " --tAccL=0"; 
}
# always full sequence length for query
$intaRNA2call .= " --qAccL=0"; 
if (defined $args{l}) { 
	$intaRNA2call .= " --tIntLenMax=".($args{l}); 
	$intaRNA2call .= " --qIntLenMax=".($args{l}); 
} else {
	$intaRNA2call .= " --tIntLenMax=0"; 
	$intaRNA2call .= " --qIntLenMax=0"; 
}
if (defined $args{o} && $args{o}==1) {
	# setup detailed v1 output
	$intaRNA2call .=" --outMode=3";
} else {
	# setup normal v1 output
	$intaRNA2call .=" --outMode=2";
}
if (defined $args{n} && $args{n}==1) {
	$intaRNA2call .= " --noSeed"
}
$intaRNA2call .= " --mode=2";
if (defined $args{u}) {
	$intaRNA2call .= " --seedMaxUP=".$args{u}
} else {
	$intaRNA2call .= " --seedMaxUP=0"
}
#$intaRNA2call .=" --seedMaxE=999"; # enable for IntaRNA v1-like seed handling
$intaRNA2call .=" --energy=F";

# call intaRNA 2
system($intaRNAbinPath."IntaRNA"." ".$intaRNA2call);
