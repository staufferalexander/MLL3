#!/usr/bin/perl                                                                                                    
#use strict;                                                                                
use warnings;
use XML::Simple;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

open (MYFILE, "/Users/staufferkm/Desktop/gencode.v19.nochr.annotation.gff3") ||die "Can't open MYFILE!";

my @gff;
while (<MYFILE>) {
    chomp;
	my @temp1 = split(/\t/, $_);
	if ($temp1[2] =~ m/gene/){
		push(@gff, $_);	
	}
    
}

print Dumper
open (OUT, ">/Users/staufferkm/Desktop/ForBedtoolsClosest_RandomGenes.gff3")||die "Can't open RandomGenes";
my $range = @gff;
#print $range;
#my @meangerps;
#my @repetitions = (1..10000);
#print Dumper \@repetitions;
#my $randmean;
my @rands;
#print $gerps[int(rand($range))];
my $numgenes = 4547; #### CHANGE THIS EVERY TIME
my $numgenes1 = $numgenes - 1;

#my $sum;
#foreach my $i(1..$#repetitions){
@rands = shuffle(@gff);
@rands1 = @rands[0..$numgenes1];
print OUT Dumper (\@rands1);


