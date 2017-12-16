#!/usr/bin/env perl

use strict;
use Getopt::Std;
use File::Basename qw( dirname );

# my $intaRNAbinPath = dirname(__FILE__)."/";
# TODO reenable above line
my $intaRNAbinPath = dirname(__FILE__)."/../src/bin/";

my %args;

if (!getopts("t:q:l:n:", \%args) or (defined $args{h} && $args{h}==1)) {
	print "CHECK PARAMTERS LIST : BUT LISTING NOT IMPLEMENTED YET";
	exit 0;
}

# TODO renenable check mandatory arguments
#if (!defined $args{t}) { die "target not specified"; };
#if (!defined $args{q}) { die "query not specified"; };
$args{t} = "aguuagucaaugaccuuuugcaccgcuuugcggugcuuuccuggaacaacaaaaugucauauacaccgaugagugaucucggacaacaaggguuguucgacaucacucggaca"; # fhlA (from figure)
$args{q} = "gaaacggagcggcaccucuuuuaacccuugaagucacugcccguuucgagaguuucucaacucgaauaacuaaagccaacgugaacuuuugcggaucuccaggauccgc"; # OxyS (from figure)
#$args{q} = "AGAATAGCAATGAACGATTATCCCTATCAAGCATTCTGACTGATAATTGCTCACAGAAACGGAGCGGCACCTCTTTTAACCCTTGAAGTCACTGCCCGTTTCGAGAGTTTCTCAACTCGAATAACTAAAGCCAACGTGAACTTTTGCGGATCTCCAGGATCCGCTTTTTTTTGCCATAAAAAA"; # OxyS (from Rfam)

$args{t} = "AUGGUUGUAAUUCCGGCAAUCUGGAGGCGUUCUGGACAGGGGUUCGAUUCCCCUCACCUCCACCA"; # accessfold example from figure
$args{q} = "GGGGGUGUACUGGUCUCGACAGGGCGGACAAAGGUGCGCAGGCAACU"; # accessfold example from figure

$args{t} = "GGCCAAGUUGUGCUCGAUUGUUCUAAGGUAACUUAGAACAGUUUGAAUGGGUUGAAUAUAGAGACCGCAUGAAUAUUC"; # CrossCatalytic paper construct
$args{q} = "GGUUCAUGUGCUCGAUUGUUACGUAAGUAACAGUUUGAAUGGGUUGAAUAUAGAGACCGCAACUUA"; # CrossCatalytic paper construct

# fill optional arguments if missing
if (!defined $args{n}) { $args{n} = 10; };
#if (!defined $args{w}) { $args{w} = 150; };
if (!defined $args{l}) { $args{l} = 100; };

# check if RNAfold accessible
if (`which RNAfold 2> /dev/null` !~ /RNAfold/) { die "RNAfold not found in PATH environment"; }

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
				." --seedBP=7 "
				." --outMode=C --outCsvCols='start1, end1, start2, end2, E, ED1, ED2, E_init, hybridDB'";

my @subopts = split /\n/, `$intaRNAbinPath/IntaRNA $intaRNAargs`;
if ( $#subopts < 2 ) { die "no alternative interactions found"; }

###########################################################################

sub max ($$) { $_[$_[0] < $_[1]] }
sub min ($$) { $_[$_[0] > $_[1]] }

###########################################################################
# for each suboptimal interaction S
###########################################################################
# get mfe interaction data for checks
my @dataMfe = split /;/, $subopts[1];
print "hdr = ".$subopts[0]."\n";
print "mfe = ".$subopts[1]."\n";
# iterate over all suboptimal interactions
for ( my $s=2; $s <= $#subopts; $s++ ) {
	# get interaction data
	my @data = split /;/, $subopts[$s];
	print "subopt = ".$subopts[$s]."\n";

	# generate IntaRNA constraints
	my $intaRNAconstr = "--tAccConstr='b:".$data[0]."-".$data[1]."' --qAccConstr='b:".$data[2]."-".$data[3]."'";

	my @constraintOut = split /\n/, `$intaRNAbinPath/IntaRNA $intaRNAargs $intaRNAconstr`;
	# check if mfe is maintained
	if ( $#subopts >= 1 ) { 
		# get constrained mfe data
		my @dataMfeConstr = split /;/, $constraintOut[1];
		print "mfe-constr = ".$constraintOut[1]."\n";
		# compute overlap
		my $overlapT = min($dataMfe[1],$dataMfeConstr[1])-max($dataMfe[0],$dataMfeConstr[0]);
		my $overlapQ = min($dataMfe[3],$dataMfeConstr[3])-max($dataMfe[2],$dataMfeConstr[2]);
		my $overlap = min($overlapT,$overlapQ) / min($dataMfe[1]-$dataMfe[0],$dataMfe[3]-$dataMfe[2]);
		# check if overlapping at all  
		if ($overlap > 0) {
			# print combined interaction details
			print "\n# compatible interactions with mfe-overlap = $overlap :\n".$subopts[$s]."\n".$constraintOut[1]."\n";
			# compute overall energy
			# generate RNAfold constraints
			my $tmp = $args{t}; $tmp =~ s/./\./g;
			my @tConstr = split //, 
			$tmp; for ( my $i=$data[0]-1; $i < $data[1]; $i++ ) { $tConstr[$i] = 'x'; }
			for ( my $i=$dataMfeConstr[0]-1; $i < $dataMfeConstr[1]; $i++ ) { $tConstr[$i] = 'x'; }
			$tmp = $args{q}; $tmp =~ s/./\./g;
			my @qConstr = split //, $tmp; 
			for ( my $i=$data[2]-1; $i < $data[3]; $i++ ) { $qConstr[$i] = 'x'; }
			for ( my $i=$dataMfeConstr[2]-1; $i < $dataMfeConstr[3]; $i++ ) { $qConstr[$i] = 'x'; }
			# compute ensemble energies of constrained ensemble
			my $targetEconstr = 0;
			$tmp = join("",@tConstr);
			if ( `echo -e "$args{t}\n$tmp" | RNAfold --maxBPspan=$args{l} --noPS -p0 -C` =~ /energy of ensemble = (\S+) kcal/ ) {	$targetEconstr = $1; }
#			print "echo -e '$args{t}\\n$tmp' | RNAfold --maxBPspan=$args{l} --noPS -p0 -C\n";
			my $queryEconstr = 0;
			$tmp = join("",@qConstr);
			if ( `echo -e "$args{q}\n$tmp" | RNAfold --maxBPspan=$args{l} --noPS -p0 -C` =~ /energy of ensemble = (\S+) kcal/ ) {	$queryEconstr = $1; }
#			print "echo -e '$args{q}\\n$tmp' | RNAfold --maxBPspan=$args{l} --noPS -p0 -C\n";
			print "combined ED1 = ".($targetEconstr-$targetEfull)." ED2 = ".($queryEconstr-$queryEfull)."\n";
			print "E = ".($data[4]-$data[5]-$data[6]
						+$dataMfeConstr[4]-$dataMfeConstr[5]-$dataMfeConstr[6]
						+$targetEconstr-$targetEfull+$queryEconstr-$queryEfull)." = E(subopt)-ED1&2(subopt)+E(mfeConstr)-ED1&2(mfeConstr)+ED1+ED2\n";
			print "E' = ".($data[4]-$data[5]-$data[6]
						+$dataMfeConstr[4]-$dataMfeConstr[5]-$dataMfeConstr[6]
						+$targetEconstr-$targetEfull+$queryEconstr-$queryEfull-$data[7])." = E(subopt)-ED1&2(subopt)+E(mfeConstr)-ED1&2(mfeConstr)+ED1+ED2-Einit\n";
			print "\n";
		}
	}
}
# generate constraint that blocks interaction site of S
# run IntaRNA
# check if mfe is maintained
# if so, 
## generate constraints blocking both interaction sites
## run RNAfold for q ant t with constraints to get subensemble energy
## compute ED and report final interaction with details


# call intaRNA 2
#system($intaRNAbinPath."IntaRNA"." ".$intaRNA2call);
