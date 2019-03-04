#!/usr/bin/env perl

############################################################################
#
# IntaRNAlink
#
#
############################################################################
# Additional dependencies:
# - RNAfold has to be accessible via PATH environment
############################################################################

use strict;
use Getopt::Std;
use File::Basename qw( dirname );

my %args;

# setup default values for input arguments
my $defCsvColsEssential = "E,start1,end1,start2,end2,ED1,ED2";
my $defCsvCols = "$defCsvColsEssential,hybridDB";
my $defL = 100;

# check if arguments ok
if (!getopts("hmt:q:l:", \%args) or (defined $args{h} && $args{h}==1)) {
	print "Available arguments:\n"
	."  -t STR\ttarget sequence (=seq1)\n"
	."  -q STR\tquery sequence (=seq2)\n"
	."  -l INT\t(opt) max. intramol. base pair span (def=$defL)\n"
	."  -m\t(opt) if given, interaction with multiple queries is modeled\n"
	."  -h\t(opt) parameter list\n"
	;
	if (defined $args{h} && $args{h}==1) {
		exit 0;
	}
	exit -1; # error in input arguments
}

# check mandatory arguments
if (!defined $args{t}) { die "ERROR: target not specified"; };
if (!defined $args{q}) { die "ERROR: query not specified"; };
# fill optional arguments if missing
if (!defined $args{l}) { $args{l} = $defL; };
my $multiQuery = "1" eq "0";
if (defined $args{m}) { $multiQuery = "1" eq "1"; };



# setup example data

# OxyS-fhlA
# IntaRNA -q AGAATAGCAATGAACGATTATCCCTATCAAGCATTCTGACTGATAATTGCTCACAGAAACGGAGCGGCACCTCTTTTAACCCTTGAAGTCACTGCCCGTTTCGAGAGTTTCTCAACTCGAATAACTAAAGCCAACGTGAACTTTTGCGGATCTCCAGGATCCGCTTTTTTTTGCCATAAAAAA -t aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca -n 2 --outmode=C --outcsvcols="id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,ED1,ED2"
#$args{t} = "aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca"; # fhlA (from figure)
#$args{q} = "gaaacggagcggcaccucuuuuaacccuugaagucacugcccguuucgagaguuucucaacucgaauaacuaaagccaacgugaacuuuugcggaucuccaggauccgc"; # OxyS (from figure)
#$args{q} = "AGAATAGCAATGAACGATTATCCCTATCAAGCATTCTGACTGATAATTGCTCACAGAAACGGAGCGGCACCTCTTTTAACCCTTGAAGTCACTGCCCGTTTCGAGAGTTTCTCAACTCGAATAACTAAAGCCAACGTGAACTTTTGCGGATCTCCAGGATCCGCTTTTTTTTGCCATAAAAAA"; # OxyS (from Rfam)

# ../src/bin/IntaRNA -q AGAATAGCAATGAACGATTATCCCTATCAAGCATTCTGACTGATAATTGCTCACAGAAACGGAGCGGCACCTCTTTTAACCCTTGAAGTCACTGCCCGTTTCGAGAGTTTCTCAACTCGAATAACTAAAGCCAACGTGAACTTTTGCGGATCTCCAGGATCCGCTTTTTTTTGCCATAAAAAA -t aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca -n 2 --outmode=C --outcsvcols="id1,start1,end1,id2,start2,end2,subseqDP,hybridDP,E,ED1,ED2" | perl IntaRNAlink.pl -t aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca -q AGAATAGCAATGAACGATTATCCCTATCAAGCATTCTGACTGATAATTGCTCACAGAAACGGAGCGGCACCTCTTTTAACCCTTGAAGTCACTGCCCGTTTCGAGAGTTTCTCAACTCGAATAACTAAAGCCAACGTGAACTTTTGCGGATCTCCAGGATCCGCTTTTTTTTGCCATAAAAAA


# accessFold example
#$args{t} = "AUGGUUGUAAUUCCGGCAAUCUGGAGGCGUUCUGGACAGGGGUUCGAUUCCCCUCACCUCCACCA"; # accessfold example from figure
#$args{q} = "GGGGGUGUACUGGUCUCGACAGGGCGGACAAAGGUGCGCAGGCAACU"; # accessfold example from figure

# CrossCatalytic example
#$args{t} = "GGCCAAGUUGUGCUCGAUUGUUCUAAGGUAACUUAGAACAGUUUGAAUGGGUUGAAUAUAGAGACCGCAUGAAUAUUC"; # CrossCatalytic paper construct
#$args{q} = "GGUUCAUGUGCUCGAUUGUUACGUAAGUAACAGUUUGAAUGGGUUGAAUAUAGAGACCGCAACUUA"; # CrossCatalytic paper construct

# Spot42-sthA
#$args{t} = "GGGATCAATTGGCTTACCCGCGATAAAATGTTACCATTCTGTTGCTTTTATGTATAAGAACAGGTAAGCCCTACCATGCCACATTCCTACGATTACGATGCCATAGTAATAGGTTCCGGCCCCGGCGGCGAAGGCGCTGCAATGGGCCTG"; # sthA (b3962) intarna-2 paper
#$args{q} = "guaggguacagagguaagauguucuaucuuucagaccuuuuacuucacguaaucggauuuggcugaauauuuuagccgccccagucaguaaugacuggggcguuuuuua"; # Spot42 intarna-2 paper




# check if RNAfold accessible
if (`which RNAfold 2> /dev/null` !~ /RNAfold/) { die "ERROR: RNAfold not found in PATH environment nor in this script's directory"; }

###########################################################################
# run RNAfold for t an q to get according ensemble energies
###########################################################################

my $Efull1 = 0;
if ( `echo $args{t} | RNAfold --maxBPspan=$args{l} --noPS -p0` =~ /energy of ensemble = (\S+) kcal/ ) {	$Efull1 = $1; }
my $Efull2 = 0;
if ( `echo $args{q} | RNAfold --maxBPspan=$args{l} --noPS -p0` =~ /energy of ensemble = (\S+) kcal/ ) {	$Efull2 = $1; }


###########################################################################
# run IntaRNA and get list of n non-overlapping suboptimal interactions
###########################################################################

# read first 2 interactions from intarna CSV output from STDIN
my @subopts = ();
while( scalar @subopts < 3 and my $line = <>) {
	$line =~ s/^\s+|\s+$//g; # trim 
	$line =~ s/\x1b\[[0-9;]*[a-zA-Z]//g; # remove color coding and escape sequences
	if ($line =~ /^\#.*/) { next; } # skip comments
	push @subopts, $line;
}

# ensure needed CSV colums are available
my $tmpString = ";".$subopts[0].";";
for my $csvCol (split /,/, $defCsvColsEssential) {
	if ( $tmpString !~ /;\s*$csvCol\s*;/) {die "ERROR: essential CSV column '$csvCol' not present in provided IntaRNA output with header '".$subopts[0]."'.\n       Ensure the following columns are present: $defCsvColsEssential"};
}

# print initial input
print $subopts[0]."\n";
if (scalar @subopts > 0) { print $subopts[1]."\n"; }
if (scalar @subopts > 1) { print $subopts[2]."\n"; }

# check if no combination possible
if (scalar @subopts < 3) {
	exit 0;
}

###########################################################################

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

# mapping of col names to col number
my %h2i;
my @dataHdr = split /;/, $subopts[0];
for ( my $i=0; $i <= $#dataHdr; $i++ ) {
	my $tmpC = $dataHdr[$i];
	$tmpC =~ s/^\s+|\s+$//g; # trim 
	$h2i{$tmpC} = $i; # store mapping
}

###########################################################################
# combine interactions
###########################################################################


# get interaction data
my @data1 = split /;/, $subopts[1];
my @data2 = split /;/, $subopts[2];

# check for overlaps
my $overlap1 = (($data1[$h2i{'start1'}] <= $data2[$h2i{'start1'}]) and ($data2[$h2i{'start1'}] <= $data1[$h2i{'end1'}]))
			or (($data1[$h2i{'start1'}] <= $data2[$h2i{'end1'}]) and ($data2[$h2i{'end1'}] <= $data1[$h2i{'end1'}]))
			or (($data2[$h2i{'start1'}] <= $data1[$h2i{'end1'}]) and ($data1[$h2i{'end1'}] <= $data2[$h2i{'end1'}]))
			or (($data2[$h2i{'start1'}] <= $data1[$h2i{'start1'}]) and ($data1[$h2i{'start1'}] <= $data2[$h2i{'end1'}]));
			
# if overlap in target, no combination possible
if ($overlap1) {
	print "DEBUG: overlap in target (seq1)\n start"." ".$data1[$h2i{'start1'}]." ".$data2[$h2i{'start1'}],$data2[$h2i{'start1'}]." ".$data1[$h2i{'end1'}];
	exit 0;
}

my $overlap2 = ($data1[$h2i{'start2'}] <= $data2[$h2i{'start2'}] and $data2[$h2i{'start2'}] <= $data1[$h2i{'end2'}])
			or ($data1[$h2i{'start2'}] <= $data2[$h2i{'end2'}] and $data2[$h2i{'end2'}] <= $data1[$h2i{'end2'}])
			or ($data2[$h2i{'start2'}] <= $data1[$h2i{'end2'}] and $data1[$h2i{'end2'}] <= $data2[$h2i{'end2'}])
			or ($data2[$h2i{'start2'}] <= $data1[$h2i{'start2'}] and $data1[$h2i{'start2'}] <= $data2[$h2i{'end2'}]);

	
########################################
# compute overall energy
########################################
# generate RNAfold constraints
my $tmp = $args{t}; $tmp =~ s/./\./g;
my @constr1 = split //, $tmp; 
for ( my $i=$data1[$h2i{'start1'}]-1; $i < $data1[$h2i{'end1'}]; $i++ ) { $constr1[$i] = 'x'; }
for ( my $i=$data2[$h2i{'start1'}]-1; $i < $data2[$h2i{'end1'}]; $i++ ) { $constr1[$i] = 'x'; }
# compute ensemble energies of constrained ensemble
my $Econstr1 = 0;
$tmp = join("",@constr1);
if ( `echo -e "$args{t}\n$tmp" | RNAfold --maxBPspan=$args{l} --noPS -p0 -C` =~ /energy of ensemble = (\S+) kcal/ ) {	$Econstr1 = $1; }
# compute ED values for combined interaction
my $ED1both = ($Econstr1-$Efull1);

# compute ED2both if not overlapping
my $ED2both = 99999;
if (not $overlap2) {
	# generate RNAfold constraints
	$tmp = $args{q}; $tmp =~ s/./\./g;
	my @constr2 = split //, $tmp; 
	for ( my $i=$data1[$h2i{'start2'}]-1; $i < $data1[$h2i{'end2'}]; $i++ ) { $constr2[$i] = 'x'; }
	for ( my $i=$data2[$h2i{'start2'}]-1; $i < $data2[$h2i{'end2'}]; $i++ ) { $constr2[$i] = 'x'; }
	# compute ensemble energies of constrained ensemble
	my $Econstr2 = 0;
	$tmp = join("",@constr2);
	if ( `echo -e "$args{q}\n$tmp" | RNAfold --maxBPspan=$args{l} --noPS -p0 -C` =~ /energy of ensemble = (\S+) kcal/ ) {	$Econstr2 = $1; }
	# compute ED values for combined interaction
	$ED2both = ($Econstr2-$Efull2);
}


if ($multiQuery) {
	########################################
	# print interaction details for multiple query molecules 
	########################################
	my @outCols = split /;/, $subopts[0];
	for (my $c=0; $c <= $#outCols; $c++) {
		#print "(".$outCols[$c].")"; 
		# corrected overall energy value
		if ($outCols[$c] eq "E") { 
			# calculate final energy of combined interaction
			my $E = # s.E - s.ED1 
					$data1[$h2i{'E'}]-$data1[$h2i{'ED1'}]
					# + c.E - c.ED1 
					+$data2[$h2i{'E'}]-$data2[$h2i{'ED1'}]
					# + ED1both
					+$ED1both;
			print $E;
		}
		elsif ($outCols[$c] eq "id1") { print $data1[$h2i{$outCols[$c]}]; }
		elsif ($outCols[$c] eq "id2") { print $data1[$h2i{$outCols[$c]}]; }
		elsif ($outCols[$c] eq "E_init") { print $data1[$h2i{$outCols[$c]}]; }
		# corrected ED1 value
		elsif ($outCols[$c] eq "ED1") { print $ED1both; }
		# default: print listed values
		else { print $data1[$h2i{$outCols[$c]}].":".$data2[$h2i{$outCols[$c]}]; }
		# separator
		if ($c < $#outCols) { print ";" };
	} 
	print "\n";
}
elsif (not $overlap2) {
	########################################
	# print interaction details for multi-site interaction of query and target 
	########################################
	my @outCols = split /;/, $subopts[0];
	for (my $c=0; $c <= $#outCols; $c++) {
		# corrected overall energy value
		if ($outCols[$c] eq "E") { 
			# calculate final energy of combined interaction
			my $E = # s.E - s.ED1 - s.ED2
					$data1[$h2i{'E'}]-$data1[$h2i{'ED1'}]-$data1[$h2i{'ED2'}]
					# + c.E - c.ED1 - c.ED2
					+$data2[$h2i{'E'}]-$data2[$h2i{'ED1'}]-$data2[$h2i{'ED2'}]
					# + ED1 + ED2
					+$ED1both +$ED2both;
			print $E;
		}
		elsif ($outCols[$c] eq "id1") { print $data1[$h2i{$outCols[$c]}]; }
		elsif ($outCols[$c] eq "id2") { print $data1[$h2i{$outCols[$c]}]; }
		elsif ($outCols[$c] eq "E_init") { print $data1[$h2i{$outCols[$c]}]; }
		# corrected ED1 value
		elsif ($outCols[$c] eq "ED1") { print $ED1both; }
		# corrected ED2 value
		elsif ($outCols[$c] eq "ED2") { print $ED2both; }
		# default: print listed values
		else { print $data1[$h2i{$outCols[$c]}].":".$data2[$h2i{$outCols[$c]}]; }
		# separator
		if ($c < $#outCols) { print ";" };
	} 
	print "\n";
}

