#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;


my @sample;
my $header = "/scratch/stauffkm/vcf/BRCA.annovar.header.txt"; #### CHANGE THIS TO A HEADER FOR YOUR FILE
my @temp;
open(AA, "$header")||die "Can't open $header!";
while(<AA>){
    chomp $_;
    @temp = split(/\t/, $_);
    push(@sample, $temp[0]);
}
close AA;

#print Dumper \@sample;
my %mut;
my $mut_file = "/scratch/stauffkm/vcf/TTNmutsinBRCA.txt"; 
#/scratch/stauffkm/vcf/KMT2AmutsinBRCA.txt
#/scratch/stauffkm/vcf/KMT2BmutsinBRCA.txt
#/scratch/stauffkm/vcf/KMT2CmutsinBRCA.txt
#/scratch/stauffkm/vcf/KMT2DmutsinBRCA.txt
#/scratch/stauffkm/vcf/PIK3AmutsinBRCA.txt
#/scratch/stauffkm/vcf/PTENmutsinBRCA.txt
#/scratch/stauffkm/vcf/KDM6AmutsinBRCA.txt
#/scratch/stauffkm/vcf/TP53mutsinBRCA.txt
#/scratch/stauffkm/vcf/TTNmutsinBRCA.txt
open(BB, "$mut_file")||die "Can't open $mut_file!";
while(<BB>){
    chomp $_;
    my @temp = split(/\t/, $_);
#   print Dumper \@temp; ### we are good up to here
    if($temp[5]=~/exonic/){
	my $pos=$temp[0].".".$temp[1]."-".$temp[2];
	$mut{$pos}{TYPE}=$temp[7];
	$mut{$pos}{GENE}=$temp[6];
	$mut{$pos}{SNV}=$temp[8];
	$mut{$pos}{AA}=$temp[9];
	$mut{$pos}{EVS6500}=$temp[10];
	$mut{$pos}{KG1}=$temp[11];
	$mut{$pos}{DB138}=$temp[12];
	$mut{$pos}{EXAC}=$temp[13];
	$mut{$pos}{SIFT}=$temp[21];
	$mut{$pos}{SIFTP}=$temp[22];
	$mut{$pos}{HVAR}=$temp[25];
	$mut{$pos}{HVARP}=$temp[26];
	$mut{$pos}{TASTE}=$temp[29];
	$mut{$pos}{TASTEP}=$temp[30];
	$mut{$pos}{MUTA}=$temp[31];
	$mut{$pos}{MUTAP}=$temp[32];
	$mut{$pos}{PROV}=$temp[35];
	$mut{$pos}{PROVP}=$temp[36];
	$mut{$pos}{GERP}=$temp[49];
	$mut{$pos}{COSMIC}=$temp[55];
	foreach my $i (68..$#temp){
		my $j=$i-2;
		my @geno = split(/:/, $temp[$i]);
		if ($geno[0] =~ /0\/1|1\/1/){
		    if ($geno[0]=~/0\/1/){
		   $mut{$pos}{$sample[$j]}{GENO}="1";
		    }
		    elsif($geno[0]=~/1\/1/){
                    $mut{$pos}{$sample[$j]}{GENO}="2";
		    }
		}
	       else{
		   $mut{$pos}{$sample[$j]}{GENO}="0";
	       }
	}
    }

#print Dumper \%mut;



my $out_file = "/scratch/stauffkm/vcf/TTNmutsinBRCAnew.txt"; ####REPLACE THIS WITH A NEW OUTFILE NAME
#/scratch/stauffkm/vcf/KMT2AmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/KMT2BmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/KMT2CmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/KMT2DmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/PIK3AmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/PTENmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/KDM6AmutsinBRCAnew.txt
#/scratch/stauffkm/vcf/TP53mutsinBRCAnew.txt
#/scratch/stauffkm/vcf/TTNmutsinBRCAnew.txt
my $out_file2 = "/scratch/stauffkm/vcf/TTNmutsinBRCAnew2";
open(OUT2, ">$out_file2")||die "Can't open $out_file2!";
open(OUT, ">$out_file")||die "Can't open $out_file!";
print OUT "pos\tgene\tSNV\tAA_change\tEVS6500\t1KG\tdb138\tExAC\tSIFT\tSIFT_pred\tPPHEN_HVAR\tPPHEN_HVAR_pred\tMTaster\tMTaster_pred\tMAssessor\tMAssessor_pred\tPROV\tPROV_pred\tGERP\tCOSMIC\t";
#print "pos\tgene\tSNV\tAA_change\tEVS6500\t1KG\tdb138\tExAC\tSIFT\tSIFT_pred\tPPHEN_HVAR\tPPHEN_HVAR_pred\tMTaster\tMTaster_pred\tMAssessor\tMAssessor_pred\tPROV\tPROV_pred\tGERP\tCOSMIC\t";
    foreach my $i (66..$#sample){
	print OUT "$sample[$i]\t";
	#print "$sample[$i]\t";
    }
    print OUT "\n";
    #print "\n";

foreach my $foo (keys %mut){

		print OUT "$foo\t$mut{$foo}{GENE}\t$mut{$foo}{SNV}\t$mut{$foo}{AA}\t$mut{$foo}{EVS6500}\t$mut{$foo}{KG1}\t$mut{$foo}{DB138}\t$mut{$foo}{EXAC}\t$mut{$foo}{SIFT}\t$mut{$foo}{SIFTP}\t$mut{$foo}{HVAR}\t$mut{$foo}{HVARP}\t$mut{$foo}{TASTE}\t$mut{$foo}{TASTEP}\t$mut{$foo}{MUTA}\t$mut{$foo}{MUTAP}\t$mut{$foo}{PROV}\t$mut{$foo}{PROVP}\t$mut{$foo}{GERP}\t$mut{$foo}{COSMIC}\t";
	#	print "$foo\t$mut{$foo}{GENE}\t$mut{$foo}{SNV}\t$mut{$foo}{AA}\t$mut{$foo}{EVS6500}\t$mut{$foo}{KG1}\t$mut{$foo}{DB138}\t$mut{$foo}{EXAC}\t$mut{$foo}{SIFT}\t$mut{$foo}{SIFTP}\t$mut{$foo}{HVAR}\t$mut{$foo}{HVARP}\t$mut{$foo}{TASTE}\t$mut{$foo}{TASTEP}\t$mut{$foo}{MUTA}\t$mut{$foo}{MUTAP}\t$mut{$foo}{PROV}\t$mut{$foo}{PROVP}\t$mut{$foo}{GERP}\t$mut{$foo}{COSMIC}\t";
		print OUT2 "$mut{$foo}{AA}\t$mut{$foo}{GERP}\t";
		foreach my $i(66..$#sample){
			print OUT "$mut{$foo}{$sample[$i]}{GENO}\t";
	#		print "$mut{$foo}{$sample[$i]}{GENO}\t";
			if ($mut{$foo}{$sample[$i]}{GENO} != "0"){
				print OUT2 "$sample[$i]\t$mut{$foo}{$sample[$i]}{GENO}\t";
			}
		}
		print OUT "\n";
		print OUT2 "\n";
	#	print "\n";
}
}
