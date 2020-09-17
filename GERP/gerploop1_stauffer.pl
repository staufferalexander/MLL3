#!/usr/bin/perl                                                                                                    
#use strict;                                                                                
use warnings;
use XML::Simple;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

open (MYFILE, "/home/stauffkm/KDM6Agerps1.txt") ||die "Can't open MYFILE!"; #### CHANGE THIS EVERYTIME

#TP53gerps1

my $mll3gerps;
while (<MYFILE>) {
    chomp;
    $mll3gerps1 .= $_;
#print $mll3gerps1;
}
my @gerps =  split(/\s/, $mll3gerps1);

open (OUT, ">/home/stauffkm/KDM6AAll_RandomGerpMeans.txt")||die "Can't open RandomGerpMeans"; #### CHANGE THIS EVERYTIME

###TP53All_RandomGerpMeans.txt

my $range = @gerps;
#print $range;
my @meangerps;
my @repetitions = (1..10000);
#print Dumper \@repetitions;
my $randmean;
my @rands;
#print $gerps[int(rand($range))];
my $numgerps = 7; #### CHANGE THIS EVERY TIME
$numgerps = $numgerps - 1;

##BRCA Mut Num with Syn SNVs (ERNEg, ERPOs, Model, All)
##TTN 68, 227, 164, 295
##PTEN 3 16 13 19
##TP53 95 94 76 189
##KMT2A 2 20 10 22
##KMT2B 3 7 5 10
##KMT2C 9 44 27 53
##KMT2D 9 17 11 26
##PIK3CA 21 168 127 189
##KDM6A 2 5 4 7 

###BRCA Mut Num without syn SNVs (ERNeg, ERPos, Model, All)
###TITIN  53, 179, 131, 232                                                                              
###pten  2, 16, 13, 18                                                                            
###tp53  92, 92, 75, 184                                                                                
###kmt2a 2, 12, 6, 14                                                                               
###kmt2b 3, 5, 3, 8                                                                               
###kmt2c 8, 38, 26, 46
###kmt2d   9, 10, 6, 19                                                                         
###pik3ca 20, 163, 124, 183
###kdm6a 2, 4, 3, 6

my $sum;
foreach my $i(1..$#repetitions){
    @rands = shuffle(@gerps);
    @rands1 = @rands[0..$numgerps];
    $sum = 0;
    foreach my $i(1 .. $#rands1) {
	$sum += $rands1[$i];
#	print "$sum\n";
    }
    $randmean = $sum/$numgerps;
    push(@meangerps, $randmean);
}
print OUT Dumper (\@meangerps);
