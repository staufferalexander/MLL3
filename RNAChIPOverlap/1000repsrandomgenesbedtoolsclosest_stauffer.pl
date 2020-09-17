#!/usr/bin/perl/
use warnings;
use strict;
use Data::Dumper;
use XML::Simple;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# open (MYFILE, "/Users/staufferkm/Desktop//Data/bedtoolsclosest/gencode.v19.nochr.annotation.gff3") ||die "Can't open MYFILE!";
# my @gff;
# while (<MYFILE>) {
#     chomp;
# 	my @temp1 = split(/\t/, $_);
# 	if ($temp1[2] =~ m/gene/){
# 		push(@gff, $_);	
# 	}
#     
# }
# print Dumper \@gff;
# foreach my $i (0..1000){
# 	my $bed = "/Users/staufferkm/Dropbox/ForBedtoolsClosest_RandomGenes_$i.gff3";
# 	open (OUT, ">$bed")||die "Can't open $bed";
# 	my $range = @gff;
# 	my $numgenes = 4547;
# 	my $numgenes1 = $numgenes - 1;
# 	my @rands = shuffle(@gff);
# 	my @rands1 = @rands[0..$numgenes1];
# 	#print OUT Dumper (\@rands1);
# 	foreach my $iter (0 .. $#rands1) {
# 		print OUT "chr$rands1[$iter]\n";	
# 	}
# }

foreach my $i (0..1000){
my $gff3 = "/Users/staufferkm/Dropbox/ForBedtoolsClosest_RandomGenes_$i.gff3";
chomp $gff3;
my @name = split(/\//, $gff3);
my @name1 = split(/\./, $name[$#name]);
my $out_dir = "/Users/staufferkm/Dropbox/RandomBedtools/";
   print "$i\n";
   my $sort = $out_dir."/".$name1[0].".sort.gff3";
   `perl /Users/staufferkm/Downloads/gff3sort-master/gff3sort.pl $gff3 > $sort`;
    `convert2bed -i gff < $gff3 > $sort`;
   print "$sort\n";
   `bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Data/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shMLL3Era.final.peaks.sort.bed -b $sort | cut -f 15 > /Users/staufferkm/Dropbox/RandomBedtools/ZME_$i.closest.txt`;
   `bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Data/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shMLL3H3K4me1.final.peaks.sort.bed -b $sort | cut -f 15   > /Users/staufferkm/Dropbox/RandomBedtools/ZMH_$i.closest.txt`;
   `bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Data/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shMLL3SP1.final.peaks.sort.bed -b $sort | cut -f 15   > /Users/staufferkm/Dropbox/RandomBedtools/ZMS_$i.closest.txt`;
   `bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Data/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shLucifEra.final.peaks.sort.bed -b $sort | cut -f 15   > /Users/staufferkm/Dropbox/RandomBedtools/ZLE_$i.closest.txt`;
   `bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Data/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shLucifH3K4me1.final.peaks.sort.bed -b $sort | cut -f 15   > /Users/staufferkm/Dropbox/RandomBedtools/ZLH_$i.closest.txt`;
   `bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Data/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shLucifSP1.final.peaks.sort.bed -b $sort | cut -f 15   > /Users/staufferkm/Dropbox/RandomBedtools/ZLS_$i.closest.txt`;
#	bedtools closest -D ref -t all -k 2 -mdb each -a /Users/staufferkm/Desktop/Sequencing/ChIP-seq/2019bedfiles/bedfiles/ZR751shMLL3Era.final.peaks.sort.bed -b /Users/staufferkm/Desktop/bedtoolsclosest/ZR751DEG_April2019.sort.gff3 > /Users/staufferkm/Desktop/DEGdistancesfromZME_bedtoolsclosest_April2019.txt
}
 
