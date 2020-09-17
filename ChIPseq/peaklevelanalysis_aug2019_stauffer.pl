#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

my %no;
my %only;
my %hash;
my %gain;
my %loss;
my %maintain;
my $pretable = "/Users/staufferkm/Desktop/Data/DEGassignedtopeaks_bycategory/PeakCategoryTable_fin_Aug2019.txt";
open(AA, "$pretable")||die "Can't open $pretable!";
while(<AA>){
    chomp $_;
    my @temp = split(/\t/, $_);
    my @temp1 = split(/\s/, $temp[1]);
    my $length = @temp1;
    my $neglength = $length*(-1);
 	my $ab = $temp[5];
 	my $dir = $temp[4];
 	my $DEG = $temp[0];
 	my $cells = $temp[6];
 	my $cat = $temp[7];
 	if ($cat =~ m/NoPeaks/){
    	if ($cells =~ m/shLucif/){
    		$hash{$DEG}{$dir}{$ab}{$cells}= [@temp1]; ## so shL = 0 0 and shM = NA
    		$loss{$DEG}{$dir}{$ab} = [@temp1]; ### so 0 0 and NA should get put here for RPL36
			push @{$no{$DEG}{$dir}{$ab}}, [@temp1];
    	}
    }
    if ($cat =~ m/OnlyPeaks/){
    	if ($cells =~ m/shMLL3/){
    		$hash{$DEG}{$dir}{$ab}{$cells}= [@temp1];
    		$gain{$DEG}{$dir}{$ab} = [@temp1];
    		push @{$only{$DEG}{$dir}{$ab}}, [@temp1];
    	}
    }
    if ($cat =~ m/LessPeaks/){
    	$hash{$DEG}{$dir}{$ab}{$cells}= [@temp1];
    }
    if ($cat =~ m/MorePeaks/){
    	$hash{$DEG}{$dir}{$ab}{$cells}= [@temp1];
    }
    if ($cat =~ m/SamePeaks/){
    	$hash{$DEG}{$dir}{$ab}{$cells}= [@temp1];
    }
}
#print Dumper \%hash;
my %niatniam;
my $shM = "shMLL3";
my $shL = "shLucif";
my %seen;
for my $deg ( sort keys %hash ){
	#print "$deg\n";
	for my $direc (keys %{$hash{$deg}}) {
		#print "$direc\n";
		for my $anti (keys %{$hash{$deg}{$direc}}) {
 			#print "$anti\n";
 			if (exists $hash{$deg}{$direc}{$anti}{$shM}){
 				unless ((exists $no{$deg}{$direc}{$anti}) or (exists $only{$deg}{$direc}{$anti})) {
					foreach my $element (@{$hash{$deg}{$direc}{$anti}{$shM}}) {
						#print "$element\n";
						my $min = $element - 300;
						my $max = $element + 300;
						if (exists $hash{$deg}{$direc}{$anti}{$shL}){
							foreach my $i (@{$hash{$deg}{$direc}{$anti}{$shL}}) {	
								if ($i >= $min && $i <= $max) { 
									unless ($seen{$deg}{$anti}{$shM}{$element} =~ m/True/) {	
										push @{$maintain{$deg}{$direc}{$anti}}, $element;
										$seen{$deg}{$anti}{$shM}{$element} = "True";
										print "$element\n";
									}
								}
							}
						}
						if ((exists $hash{$deg}{$direc}{$anti}{$shL}) && (!($seen{$deg}{$anti}{$shM}{$element} =~ m/True/))) {   
							push @{$gain{$deg}{$direc}{$anti}}, $element;
						}	
	 				}
				}		
			}			
			if (exists $hash{$deg}{$direc}{$anti}{$shL}) {
				unless ((exists $no{$deg}{$direc}{$anti}) or (exists $only{$deg}{$direc}{$anti})) {
					foreach my $el (@{$hash{$deg}{$direc}{$anti}{$shL}}) {
						my $mini = $el - 300;
						my $maxi = $el + 300;
						if (exists $hash{$deg}{$direc}{$anti}{$shM}) {
							foreach my $e (@{$hash{$deg}{$direc}{$anti}{$shM}}) {	
								if ($e >= $mini && $e <= $maxi) { 
									unless ($seen{$deg}{$anti}{$shL}{$el} =~ m/True/) {	
										push @{$maintain{$deg}{$direc}{$anti}}, $el;
										$seen{$deg}{$anti}{$shL}{$el} = "True";
									}
								}
							}
						}		
						if ((exists $hash{$deg}{$direc}{$anti}{$shM}) && (!($seen{$deg}{$anti}{$shL}{$el} =~ m/True/))) {   
								push @{$loss{$deg}{$direc}{$anti}}, $el;
						}
					}
				}	
			}
		}
	}
}

#print Dumper \%seen;
#print Dumper \%loss;
#print Dumper \%maintain;
#print Dumper \%gain;
my $peaktable = "/Users/staufferkm/Desktop/PeakLevelAnalysisTable1_Aug2019.txt";
my $peaktable = "/Users/staufferkm/Desktop/PeakLevelAnalysisTable_meltfixed_Aug2019.txt";
open (OUT, ">$peaktable") || die "Can't open $peaktable!";
my $H3K4me1 = "H3K4me1";
my $SP1 = "SP1";
my $ER = "ERa";
my $up = "Up";
my $down = "Down";
## need to structure this so that we go through hash i think . . . so it prints out DEG up, maintained ERa H S, Lost ERA H S, Gain ERa H S; DEG down 
for my $deg ( sort keys %hash ){
	if ($hash{$deg}{$up}) {
		if (exists $maintain{$deg}{$up}{$ER}) {
			#print OUT "$deg\tUp\tERa maintain:  ";
			foreach my $k (@{$maintain{$deg}{$up}{$ER}}) {
				print OUT "$deg\tUp\tERa maintain\t$k\n";
			}
		#print OUT "\n";
		}
		if (exists $loss{$deg}{$up}{$ER}) {
			#print OUT "$deg\tUp\tERa loss:  ";
			foreach my $j (@{$loss{$deg}{$up}{$ER}}) {	
				if (exists $maintain{$deg}{$up}{$ER}[$j]) {
					print OUT "$deg\tUp\tERa loss\t$j\tDOUBLEDIPPED\n";  #*** only zeros are showing up as double dipped
				}
				if (!(exists $maintain{$deg}{$up}{$ER}[$j])) {
					print OUT "$deg\tUp\tERa loss\t$j\n";
				}
			}
			#print OUT "\n";
		}
		if (exists $gain{$deg}{$up}{$ER}) {
			#print OUT "$deg\tUp\tERa gain:  ";
			foreach my $g (@{$gain{$deg}{$up}{$ER}}) {	
				if (exists $maintain{$deg}{$up}{$ER}[$g]) {
					print OUT "$deg\tUp\tERa gain\t$g\tDOUBLEDIPPED\n";
				}
				if (!(exists $maintain{$deg}{$up}{$ER}[$g])) {
					print OUT "$deg\tUp\tERa gain\t$g\n";
				}
			}
			#print OUT "\n";
		}
	}
	if ($hash{$deg}{$down}) {
		if (exists $maintain{$deg}{$down}{$ER}) {
			#print OUT "$deg\tDown\tERa maintain:  ";
			foreach my $x (@{$maintain{$deg}{$down}{$ER}}) {
				print OUT "$deg\tDown\tERa maintain\t$x\n";
			}
		#print OUT "\n";
		}
		if (exists $loss{$deg}{$down}{$ER}) {
			#print OUT "$deg\tDown\tERa loss:  ";
			foreach my $y (@{$loss{$deg}{$down}{$ER}}) {	
				if (exists $maintain{$deg}{$down}{$ER}[$y]){
					print OUT "$deg\tDown\tERa loss\t$y\tDOUBLEDIPPED\n";
				}
				if (!(exists $maintain{$deg}{$down}{$ER}[$y])){
					print OUT "$deg\tDown\tERa loss\t$y\n";
				}
			}
			#print OUT "\n";
		}
		if (exists $gain{$deg}{$down}{$ER}) {
			#print OUT "$deg\tDown\tERa gain:  ";
			foreach my $z (@{$gain{$deg}{$down}{$ER}}) {	
				if (exists $maintain{$deg}{$down}{$ER}[$z]){
					print OUT "$deg\tDown\tERa gain\t$z\tDOUBLEDIPPED\n";
				}
				if (!(exists $maintain{$deg}{$down}{$ER}[$z])){
					print OUT "$deg\tDown\tERa gain\t$z\n";
				}
			}
			#print OUT "\n";
		}
	}
}	

###need to ask the following:
###		distance to promoter (the number of the peak assigned to the DEG). For Upreg genes --> for ERa --> for lost peaks
###																									   --> for gained peaks
###																									   --> for maintained peaks
###																		 For downreg genes --> for ERa --> for lost peaks
###																									   --> for gained peaks
###																									   --> for maintained peaks
###    motif analysis  	(basically need to just print out the genomic locations for each gene for each category, so I can then make a file for MEME?)		
###																		   For Upreg genes --> for ERa --> for lost peaks
###																									   --> for gained peaks
###																									   --> for maintained peaks
###																		 For downreg genes --> for ERa --> for lost peaks
###																									   --> for gained peaks
###																									   --> for maintained peaks
###    H3K4me1 & SP1 association

### this will not print out duplicates (if shMLL3 has six 0s for a gene, and shLucif has two 0s, it will only print out two maintain 0s)