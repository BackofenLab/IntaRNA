#!/usr/bin/env perl

############################################################################
#
# IntaRNA-fuse
#
# Computes non-overlapping suboptimal interaction sites and checks whether
# they are compatible withe the mfe interaction. That is whether a constrained
# IntaRNA prediction (where the suboptimal site is marked as blocked) shows
# an overlap of its constrained mfe with the unconstrained mfe interaction.
# If so, the correct overall energy of the fused interaction is computed and
# a joint interaction CSV information is reported.
#
############################################################################
# Additional dependencies:
# - RNAfold has to be accessible via PATH environment
############################################################################

use strict;
use Getopt::Std;
use File::Basename qw( dirname );

my $intaRNAbinPath = dirname(__FILE__)."/";
if ( `$intaRNAbinPath/IntaRNA -h 2>/dev/null` !~ /RNA-RNA/ ) {
	# check if in compilation environment and not installed yet
	if ( `$intaRNAbinPath/../src/bin/IntaRNA -h 2>/dev/null` =~ /RNA-RNA/ ) {
		$intaRNAbinPath = dirname(__FILE__)."/../src/bin/";
	} else {
		die "ERROR: IntaRNA not within this script's directory!";
	}
}

my %args;

# setup default values for input arguments
my $defN = 10;
my $defL = 100;
my $defCsvColsEssential = "E,start1,end1,start2,end2,ED1,ED2,E_init";
my $defCsvCols = "$defCsvColsEssential,hybridDB";

# check if arguments ok
if (!getopts("ht:q:l:n:o:", \%args) or (defined $args{h} && $args{h}==1)) {
	print "Available arguments:\n"
	."  -t\ttarget sequence (=seq1)\n"
	."  -q\tquery sequence (=seq2)\n"
	."  -n\t(opt) max. number of subopts (def=$defN)\n"
	."  -l\t(opt) max. intramol. base pair span (def=$defL)\n"
	."  -o\t(opt) CSV columns to report (def=$defCsvCols)\n"
	."  -h\t(opt) parameter list\n"
	;
	if (defined $args{h} && $args{h}==1) {
		exit 0;
	}
	exit -1; # error in input arguments
}

# setup example data

# OxyS-fhlA
$args{t} = "aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca"; # fhlA (from figure doi:10.1006/jmbi.2000.3942:Fig7)
$args{q} = "gaaacggagcggcaccucuuuuaacccuugaagucacugcccguuucgagaguuucucaacucgaauaacuaaagccaacgugaacuuuugcggaucuccaggauccgc"; # OxyS (from figure doi:10.1006/jmbi.2000.3942:Fig7) 
$args{q} = "AGAAUAGCAAUGAACGAUUAUCCCUAUCAAGCAUUCUGACUGAUAAUUGCUCACAGAAACGGAGCGGCACCUCUUUUAACCCUUGAAGUCACUGCCCGUUUCGAGAGUUUCUCAACUCGAAUAACUAAAGCCAACGUGAACUUUUGCGGAUCUCCAGGAUCCGCUUUUUUUUGCCAUAAAAAA"; # OxyS (from Rfam)

# Spot42 example https://www.embopress.org/cgi/doi/10.1038/emboj.2012.52

# accessFold example
#$args{t} = "AUGGUUGUAAUUCCGGCAAUCUGGAGGCGUUCUGGACAGGGGUUCGAUUCCCCUCACCUCCACCA"; # accessfold example from figure
#$args{q} = "GGGGGUGUACUGGUCUCGACAGGGCGGACAAAGGUGCGCAGGCAACU"; # accessfold example from figure

# CrossCatalytic example
#$args{t} = "GGCCAAGUUGUGCUCGAUUGUUCUAAGGUAACUUAGAACAGUUUGAAUGGGUUGAAUAUAGAGACCGCAUGAAUAUUC"; # CrossCatalytic paper construct
#$args{q} = "GGUUCAUGUGCUCGAUUGUUACGUAAGUAACAGUUUGAAUGGGUUGAAUAUAGAGACCGCAACUUA"; # CrossCatalytic paper construct

# Spot42-sthA
#$args{t} = "GGGATCAATTGGCTTACCCGCGATAAAATGTTACCATTCTGTTGCTTTTATGTATAAGAACAGGTAAGCCCTACCATGCCACATTCCTACGATTACGATGCCATAGTAATAGGTTCCGGCCCCGGCGGCGAAGGCGCTGCAATGGGCCTG"; # sthA (b3962) intarna-2 paper
#$args{q} = "guaggguacagagguaagauguucuaucuuucagaccuuuuacuucacguaaucggauuuggcugaauauuuuagccgccccagucaguaaugacuggggcguuuuuua"; # Spot42 intarna-2 paper

# check mandatory arguments
if (!defined $args{t}) { die "ERROR: target not specified"; };
if (!defined $args{q}) { die "ERROR: query not specified"; };
# fill optional arguments if missing
if (!defined $args{n}) { $args{n} = $defN; };
if (!defined $args{l}) { $args{l} = $defL; };
if (!defined $args{o}) { $args{o} = $defCsvCols; };

# ensure needed CSV colums are available
for my $csvCol (split /,/, $defCsvColsEssential) {
	my $tmpString = ",".$args{o}.",";
	if ( $tmpString !~ /,\s*$csvCol\s*,/) {die "ERROR: essential CSV column '$csvCol' not present in argument '-c'"};
}


# check if RNAfold accessible
if (`which RNAfold 2> /dev/null` !~ /RNAfold/) { die "ERROR: RNAfold not found in PATH environment nor in this script's directory"; }

###########################################################################
# run RNAfold for t an q to get according ensemble energies
###########################################################################

my $targetEfull = 0;
if ( `echo $args{t} | RNAfold --maxBPspan=$args{l} --noPS -p0` =~ /energy of ensemble = (\S+) kcal/ ) {	$targetEfull = $1; }
my $queryEfull = 0;
if ( `echo $args{q} | RNAfold --maxBPspan=$args{l} --noPS -p0` =~ /energy of ensemble = (\S+) kcal/ ) {	$queryEfull = $1; }


###########################################################################
# run IntaRNA and get list of n non-overlapping suboptimal interactions
###########################################################################

my $intaRNAargs = "-t $args{t} --tAccW=0 --tAccL=$args{l}"
				." -q $args{q} --qAccW=0 --qAccL=$args{l}"
				." --outOverlap=N -n $args{n} "
				." --outMode=C --outCsvCols='$args{o}'";

my @subopts = split /\n/, `$intaRNAbinPath/IntaRNA $intaRNAargs`;

###########################################################################

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

###########################################################################
# for each suboptimal interaction S
###########################################################################
# get mfe interaction data for checks
my @dataMfe = split /;/, $subopts[1];
print $subopts[0]."\n";
# mapping of col names to col number
my %h2i;
my @dataHdr = split /;/, $subopts[0];
for ( my $i=0; $i <= $#dataHdr; $i++ ) {
	my $tmpC = $dataHdr[$i];
	$tmpC =~ s/^\s+|\s+$//g; # trim 
	$h2i{$tmpC} = $i; # store mapping
}
print $subopts[1]."\n";
# iterate over all suboptimal interactions
for ( my $s=2; $s <= $#subopts; $s++ ) {
	# get interaction data
	my @data = split /;/, $subopts[$s];
	print $subopts[$s]."\n";

	# generate IntaRNA constraints
	my $intaRNAconstr = "--tAccConstr='b:".$data[$h2i{'start1'}]."-".$data[$h2i{'end1'}]."' --qAccConstr='b:".$data[$h2i{'start2'}]."-".$data[$h2i{'end2'}]."'";

	my @constraintOut = split /\n/, `$intaRNAbinPath/IntaRNA $intaRNAargs $intaRNAconstr`;
	# check if mfe is maintained
	if ( $#subopts >= 1 ) { 
		# get constrained mfe data
		my @dataMfeConstr = split /;/, $constraintOut[1];
#		print "mfe-constr = ".$constraintOut[1]."\n";
		# compute overlap
		my $overlapT = min($dataMfe[$h2i{'end1'}],$dataMfeConstr[$h2i{'end1'}])-max($dataMfe[$h2i{'start1'}],$dataMfeConstr[$h2i{'start1'}]);
		my $overlapQ = min($dataMfe[$h2i{'end2'}],$dataMfeConstr[$h2i{'end2'}])-max($dataMfe[$h2i{'start2'}],$dataMfeConstr[$h2i{'start2'}]);
		my $overlap = min($overlapT,$overlapQ) / min($dataMfe[$h2i{'end1'}]-$dataMfe[$h2i{'start1'}],$dataMfe[$h2i{'end2'}]-$dataMfe[$h2i{'start2'}]);
		# check if overlapping at all  
		if ($overlapT > 0 and $overlapQ > 0) {
#		if ($overlap > 0) {
#			print "\n# compatible interactions with mfe-overlap = $overlap :\n".$subopts[$s]."\n".$constraintOut[1]."\n";
			########################################
			# compute overall energy
			########################################
			# generate RNAfold constraints
			my $tmp = $args{t}; $tmp =~ s/./\./g;
			my @tConstr = split //, 
			$tmp; for ( my $i=$data[$h2i{'start1'}]-1; $i < $data[$h2i{'end1'}]; $i++ ) { $tConstr[$i] = 'x'; }
			for ( my $i=$dataMfeConstr[$h2i{'start1'}]-1; $i < $dataMfeConstr[$h2i{'end1'}]; $i++ ) { $tConstr[$i] = 'x'; }
			$tmp = $args{q}; $tmp =~ s/./\./g;
			my @qConstr = split //, $tmp; 
			for ( my $i=$data[$h2i{'start2'}]-1; $i < $data[$h2i{'end2'}]; $i++ ) { $qConstr[$i] = 'x'; }
			for ( my $i=$dataMfeConstr[$h2i{'start2'}]-1; $i < $dataMfeConstr[$h2i{'end2'}]; $i++ ) { $qConstr[$i] = 'x'; }
			# compute ensemble energies of constrained ensemble
			my $targetEconstr = 0;
			$tmp = join("",@tConstr);
			if ( `echo -e "$args{t}\n$tmp" | RNAfold --maxBPspan=$args{l} --noPS -p0 -C` =~ /energy of ensemble = (\S+) kcal/ ) {	$targetEconstr = $1; }
			my $queryEconstr = 0;
			$tmp = join("",@qConstr);
			if ( `echo -e "$args{q}\n$tmp" | RNAfold --maxBPspan=$args{l} --noPS -p0 -C` =~ /energy of ensemble = (\S+) kcal/ ) {	$queryEconstr = $1; }
			# compute ED values for combined interaction
			my $ED1 = ($targetEconstr-$targetEfull);
			my $ED2 = ($queryEconstr-$queryEfull);
			# calculate final energy of combined interaction
			my $E = # s.E - s.ED1 - s.ED2
					$data[$h2i{'E'}]-$data[$h2i{'ED1'}]-$data[$h2i{'ED2'}]
					# + c.E - c.ED1 - c.ED2
					+$dataMfeConstr[$h2i{'E'}]-$dataMfeConstr[$h2i{'ED1'}]-$dataMfeConstr[$h2i{'ED2'}]
					# + ED1 + ED2
					+$ED1 +$ED2;
			# calculate energy minus one E_init term
#			my $Eprime = $E - $data[$h2i{'E_init'}]);

			########################################
			# print combined interaction details 
			########################################
			my @outCols = split /,/, $args{o};
			for (my $c=0; $c <= $#outCols; $c++) {
				if ($outCols[$c] eq "E") { print $E; }
				elsif ($outCols[$c] eq "ED1") { print $ED1; }
				elsif ($outCols[$c] eq "ED2") { print $ED2; }
				else { print $data[$h2i{$outCols[$c]}].":".$dataMfeConstr[$h2i{$outCols[$c]}]; }
				if ($c < $#outCols) { print ";" };
			} 
			print "\n";
		}
	}
}

