#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Math::Round;

### first section
# my %ogs;
# my $convert = "/Users/staufferkm/Desktop/DEGassocwithERapeaks_Aug2019.txt"; 
# open(AA, "$convert")||die "Can't open $convert!";
# while(<AA>){
#     chomp $_;
#     #my @temp = split(/\t/, $_);
# 	#push(@ogs, $temp[0]);
# 	$ogs{$_}++;
#     }
# #print Dumper \%ogs;
# 
# my %melt;
# my $dir;
# my $cat;
# my $pretable = "/Users/staufferkm/Desktop/PeakLevelAnalysisTable_melt_Aug2019.txt";
# open(AA, "$pretable")||die "Can't open $pretable!";
# while(<AA>){
#     chomp $_;
#     my @temp = split(/\t/, $_);
#  	my $DEG = $temp[0];
#  	$dir = $temp[1];
#  	$cat = $temp[2];
#  	my $peak = $temp[3];
#     push @{$melt{$DEG}{$dir}{$cat}}, $peak;
# 	}
# #print Dumper \%melt;
# 
# my %seen;
# my $NA = "NA";
# my $out_table = "/Users/staufferkm/Desktop/DEGassocwERa_peaksforMEME_Aug2019.txt"; 
# open(OUT, ">$out_table")||die "Can't open $out_table!";
# my $gene_table = "/Users/staufferkm/Desktop/bedtoolsclosest/gencode.v19.nochr.annotation.gff3";
# open(BB, "$gene_table")||die "Can't open $gene_table!";
# while(<BB>){
#     chomp $_;
#     my @temp1 = split(/\t/, $_);
#     my $genebeg = $temp1[3];
#     my $geneend = $temp1[4];
# 	my @temp2 = split(/\;/, $temp1[8]);
# 	my @temp3 = split(/\=/, $temp2[5]);
# 	my $gene = $temp3[1];
# 	if ($ogs{$gene}){
# 		my $tss = $temp1[3];
# 		my $chr = $temp1[0];
# 		for my $direc ( keys %{$melt{$gene}} ){
# 			for my $categ (keys %{$melt{$gene}{$direc}}) {
# 				foreach my $i (@{$melt{$gene}{$direc}{$categ}}) {
# 					unless ($melt{$gene}{$direc}{$categ}[$i] =~ $NA) {
# 						my $begin = $i - 150;
# 						my $end = $i + 150;
# 						my $b = $tss + $begin;
# 						my $e = $tss + $end;
# 						if ($seen{$gene}{$direc}{$categ}{$i}) {		
# 							print OUT "$chr\t$genebeg\t$geneend\t$b\t$e\t$gene\t$direc\t$categ\t$i\n";	
# 							#print "$chr\t$gene\t$categ\t$i\n";
# 						}
#  						if (!($seen{$gene}{$direc}{$categ}{$i})) {	
#  							$seen{$gene}{$direc}{$categ}{$i} = "True";		
#  							print OUT "$chr\t$genebeg\t$geneend\t$b\t$e\t$gene\t$direc\t$categ\t$i\n";
#  						}	
#  					}
#  				}
#  			}	
#  		}
# 	}
# }

### second section
### checking for gene duplicates at different positions in the genome in the output file
# my %check;
# my $gene;
# my $checkfile = "/Users/staufferkm/Desktop/DEGassocwERa_peaksforMEME_Aug2019.txt"; 
# open(AA, "$checkfile")||die "Can't open $checkfile!";
# while(<AA>){
#     chomp $_;
# 	my @temp = split(/\t/, $_);
# 	$gene = $temp[5];
# 	my $loc = "$temp[0]\t$temp[1]\t$temp[2]";
# 	if (!($loc ~~ @{$check{$gene}})){
# 		push @{$check{$gene}}, $loc;
# 	}
# }
# #print Dumper \%check;
# for my $i (keys %check) {
# 	my $length = @{$check{$i}};
# 	#print "$length\n";
# 	if ($length > 1) {
# 		foreach my $k (@{$check{$i}}){
# 			print "$i\t$k\n";
# 		}
# 	}
# }

###third section -## update - I have decided to do this a different way! see fourth section below
#### need to somehow match those up to /Users/staufferkm/Desktop/DEGassocwERa_peaksforMEME_Aug2019.txt
### the list of duplicate genes at diff genomic locations is in /Users/staufferkm/Desktop/DuplicateGeneNames_peaksforMEME_Aug2019.txt
# my %dups;
# my $dup = "/Users/staufferkm/Desktop/DuplicateGeneNames_namesonly_Aug2019.txt"; 
# open(AA, "$dup")||die "Can't open $dup!";
# while(<AA>){
# 	chomp $_;
# 	$dups{$_}++;
# }
# #print Dumper \%dups;
# 
# my %table;
# my $pretable = "/Users/staufferkm/Desktop/DEGassocwERa_peaksforMEME_Aug2019.txt";
# open(BB, "$pretable")||die "Can't open $pretable!";
# while(<BB>){
# 	chomp $_;
#     my @temp = split(/\t/, $_);
#   	my $chrom = $temp[0];
#   	my $deg = $temp[5];
#   	my $pos = $temp[8];
#   	my $GS = $temp[1];
#   	my $GE = $temp[2];
#   	my $G_id = "$chrom\_"."$GS-"."$GE";
#     $table{$deg}{$pos}{$G_id} = $_;
# }
# #print Dumper \%table;
# 
# my %orig;
# my $gene;
# my $peak;
# my $gene_id;
# my $orig = "/Users/staufferkm/Desktop/AssignedZE_1MilfromDEG.txt";
# open(CC, "$orig")||die "Can't open $orig!";
# while(<CC>){
# 	chomp $_;
# 	my @temp1 = split(/\t/, $_);
# 	my $t12 = $temp1[12];
# 	my @temp2 = split(/\;/, $t12);
# 	my $t5 = $temp2[5];
#  	my @temp3 = split(/\=/, $t5);
#  	$gene = $temp3[1];
#  	$peak = $temp1[13];
#  	my $genes = $temp1[7];
#  	my $genee = $temp1[8];
#  	my @temp4 = split('r', $temp1[0]);
#  	my $chr = $temp4[1];
#  	$gene_id = "$chr\_"."$genes-"."$genee";
#  	$orig{$gene}{$peak} = $gene_id;
# }
# #print Dumper \%orig;
# 	
# my $out_table = "/Users/staufferkm/Desktop/DEGassocwERa_peaksforMEME_nodups_Sept2019.txt"; 
# open(OUT, ">$out_table")||die "Can't open $out_table!";
# for my $g (keys %table){
# 	if(exists($dups{$g})) {
# 		for my $i (keys %{$orig{$g}}){
# 			#print OUT "$orig{$g}{$i}\tDUP\n";
# 			print OUT "$table{$g}{$i}{$orig{$g}{$i}}\tDUP\n";
# 		}
# 	}
# 	if(!(exists($dups{$g}))) {
# 		for my $p (keys %{$table{$g}}) {
# 			for my $t (keys %{$table{$g}{$p}}) {
# 				print OUT "$table{$g}{$p}{$t}\tORIG\n";
# 			}
# 		}
# 	}
# }

####fourth section - I've got the correct number of peaks now (diff from the AssignedZE file because AGAP4 and some other A gene were ncRNAs that aren't included)
###so now I just need to create the bedfile based off of the genomic locations of the DEG
# my %melt;
# my $dir;
# my $cat;
# my $DEG;
# my $pretable = "/Users/staufferkm/Desktop/PeakLevelAnalysisTable_meltfixed_Aug2019.txt";
# open(AA, "$pretable")||die "Can't open $pretable!";
# while(<AA>){
#     chomp $_;
#     my @temp = split(/\t/, $_);
#  	$DEG = $temp[0];
#  	$dir = $temp[1];
#  	$cat = $temp[2];
#  	my $peak = $temp[3];
#     push @{$melt{$DEG}{$dir}{$cat}}, $peak;
# 	}
# print Dumper \%melt;
# 
# my %orig;
# my $gene;
# my $peak;
# my $peak_id;
# my $orig = "/Users/staufferkm/Desktop/FilesToSortLater/AssignedZE_1MilfromDEG.txt";
# open(CC, "$orig")||die "Can't open $orig!";
# while(<CC>){
# 	chomp $_;
# 	my @temp1 = split(/\t/, $_);
# 	my $t12 = $temp1[12];
# 	print "$t12\n";
# 	my @temp2 = split(/\;/, $t12);
# 	my $t5 = $temp2[5];
#  	print "$t5\n";
#  	my @temp3 = split(/\=/, $t5);
#  	my $gene = $temp3[1];
#  	my $peak = $temp1[13];
#  	my $peaks = $temp1[1];
#  	my $peake = $temp1[2];
#  	my @temp4 = split('r', $temp1[0]);
#  	my $chr = $temp4[1];
#  	print "$chr\n";
#  	my $peak_id = "$chr\_"."$peaks\_"."$peake";
#  	$orig{$gene}{$peak} = $peak_id;
# }
# print Dumper \%orig;
# 
# my $out_table = "/Users/staufferkm/Desktop/DEGassocwERa_peaksforMEME_final_Sept2019.bed"; 
# open(OUT, ">$out_table")||die "Can't open $out_table!";
# for my $k (keys %melt){
# 	for my $l (keys %{$melt{$k}}) {
# 		for my $m (keys %{$melt{$k}{$l}}) {
# 			foreach my $n (@{$melt{$k}{$l}{$m}}) {
# 				need to make a 500bp bed file here
# 				my @temp5 = split(/\_/, $orig{$k}{$n});
# 				my $chrom = $temp5[0];
# 				my $peakstart = $temp5[1];
# 				my $peakend = $temp5[2];
# 				my $middle = ((($peakend - $peakstart)/2) + $peakstart);
# 				my $start = $middle - 250;
# 				my $end = $middle + 250;
# 				print OUT "$chrom\t$peakstart\t$peakend\t$k\t$l\t$m\n";
# 			}
# 		}
# 	}
# }

# check slack and email for perl script and any command line code
# 
# if within the gene body --> 0 
# if peak after gene body --> gene end - peak start 
# if peak before gene body --> gene start - peak end	

###extra add on because I need to make all the peaks on the background bed file 500bp long as well so that they match up with the input files
# my $pretable = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ZR751shLshMERa_backgroundforGREAT.bed";
# my $out_table = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ZR751shLshMERa_backgroundforGREAT_500bp.bed"; 
# open(OUT, ">$out_table")||die "Can't open $out_table!";
# open(AA, "$pretable")||die "Can't open $pretable!";
# while(<AA>){
# 	my @temp = split(/\t/, $_);
#  	my $beg = $temp[1];
#  	my $end = $temp[2];
#  	my $mid = ((round($end - $beg)/2) + $beg);
# 	my $start = $mid - 250;
# 	my $last = $mid + 250;
# 	print OUT "$temp[0]\t$start\t$last\t$temp[3]\n";
# }

##ok now the foreground is not a subset of the background on GREAT . . . so need to adjust the background to match the foreground
my $down = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ERapeaks_assocwithDownDEG_forMEME_Sept2019.clip.bed";
my $up = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ERapeaks_assocwithUpDEG_forMEME_Sept2019.clip.bed";
my $lost = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ERapeakslost_assocwithDEG_forMEME_Sept2019.clip.bed";
my $gain = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ERapeaksgained_assocwithDEG_forMEME_Sept2019.clip.bed";
my $main = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ERapeaksmaintained_assocwithDEG_forMEME_Sept2019.clip.bed";
my $bg = "/Users/staufferkm/Desktop/PeakAssocDEG_MEME_Sept2019/ZR751shLshMERa_backgroundforGREAT_500bp_fixed.bed";

my %hdown;
open(AA, "$down")||die "Can't open $down!";
while(<AA>){
	chomp $_;
	my @temp = split(/\t/, $_);
	my $chr = $temp[0];
	push @{$hdown{$chr}}, $_;
}

my %hup;
open(BB, "$up")||die "Can't open $up!";
while(<BB>){
	chomp $_;
	my @temp1 = split(/\t/, $_);
	my $chr1 = $temp1[0];
	push @{$hup{$chr1}}, $_;
}

my %hlost;
open(CC, "$lost")||die "Can't open $lost!";
while(<CC>){
	chomp $_;
	my @temp2 = split(/\t/, $_);
	my $chr2 = $temp2[0];
	push @{$hlost{$chr2}}, $_;
}

my %hgain;
open(DD, "$gain")||die "Can't open $gain!";
while(<DD>){
	chomp $_;
	my @temp3 = split(/\t/, $_);
	my $chr3 = $temp3[0];
	push @{$hgain{$chr3}}, $_;
}

my %hmain;
open(EE, "$main")||die "Can't open $main!";
while(<EE>){
	chomp $_;
	my @temp4 = split(/\t/, $_);
	my $chr4 = $temp4[0];
	push @{$hmain{$chr4}}, $_;
}

my %hbg;
open(FF, "$bg")||die "Can't open $bg!";
while(<FF>){
	chomp $_;
	my @temp5 = split(/\t/, $_);
	my $chr5 = $temp5[0];
	push @{$hbg{$chr5}}, $_;
}
#print Dumper \%hbg;
#print Dumper \%hmain;

my %missing;
my $out_table = "/Users/staufferkm/Desktop/missingpeaksinbackground_try2_Sept2019.txt"; 
open(OUT, ">$out_table")||die "Can't open $out_table!";
for my $i (keys %hdown) {
	foreach my $j (@{$hdown{$i}}) {
		if ($hbg{$i}[$j]){
			print "Here\n";
		}
		else {
			$missing{$j}++;
		}
	}
}

for my $k (keys %hup) {
	foreach my $l (@{$hup{$k}}) {
		if ($hbg{$k}[$l]){
			print "Here\n";
		}
		else {
			$missing{$l}++;
		}
	}
}

for my $m (keys %hlost) {
	foreach my $n (@{$hlost{$m}}) {
		if (!($hbg{$m}[$n])){
			print "Here\n";
		}
		else {
			$missing{$n}++;
		}
	}
}

for my $o (keys %hgain) {
	foreach my $p (@{$hgain{$o}}) {
		if (!($hbg{$o}[$p])){
			print "Here\n";
		}
		else {
			$missing{$p}++;
		}
	}
}

for my $q (keys %hmain) {
	foreach my $r (@{$hmain{$q}}) {
		if (!($hbg{$q}[$r])){
			print "Here\n";
		}
		else {
			$missing{$r}++;
		}
	}
}
#print Dumper \%missing;

for my $s (keys %missing) {
	print OUT "$s\n";
}