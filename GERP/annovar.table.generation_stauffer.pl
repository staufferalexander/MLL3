#!/usr/bin/perl
#use strict;
#use warnings;
use Data::Dumper;


my @sample;
my $header = "/data/strickt2/3102/final.annovar.header.txt";
open(AA, "$header")||die "Can't open $header!";
while(<AA>){
    chomp $_;
    my @temp = split(/\t/, $_);
    push(@sample, $temp[0]);
}
close AA;

print Dumper \@sample;

my %mut;
my $mut_file = "/data/strickt2/3102/3102.GCVF.indel.snpeff.cosmic.annovar.102115.haplo.genotypes.hg19_multianno.txt";
open(BB, "$mut_file")||die "Can't opne $mut_file!";
while(<BB>){
    chomp $_;
    my @temp = split(/\t/, $_);
    if($temp[5]=~/exonic/){
	foreach my $i (67..$#temp){
	    my $j=$i-2;
	    #print "$j\t$sample[$j]\n";
	    my @geno = split(/:/, $temp[$i]);
	    if ($geno[0] =~ /0\/1|1\/1/){
		my $pos=$temp[0].".".$temp[1]."-".$temp[2];
		#$mut{$sample[$i]}{$pos}{TYPE}=$temp[7];
		$mut{$sample[$j]}{$pos}{GENE}=$temp[6];
		$mut{$sample[$j]}{$pos}{SNV}=$temp[8];
		$mut{$sample[$j]}{$pos}{AA}=$temp[9];
		$mut{$sample[$j]}{$pos}{EVS6500}=$temp[10];
		$mut{$sample[$j]}{$pos}{KG1}=$temp[11];
		$mut{$sample[$j]}{$pos}{DB138}=$temp[12];
		$mut{$sample[$j]}{$pos}{EXAC}=$temp[13];
		$mut{$sample[$j]}{$pos}{SIFT}=$temp[21];
		$mut{$sample[$j]}{$pos}{SIFTP}=$temp[22];
		$mut{$sample[$j]}{$pos}{HVAR}=$temp[25];
		$mut{$sample[$j]}{$pos}{HVARP}=$temp[26];
		$mut{$sample[$j]}{$pos}{TASTE}=$temp[29];
		$mut{$sample[$j]}{$pos}{TASTEP}=$temp[30];
		$mut{$sample[$j]}{$pos}{MUTA}=$temp[31];
		$mut{$sample[$j]}{$pos}{MUTAP}=$temp[32];
		$mut{$sample[$j]}{$pos}{PROV}=$temp[35];
		$mut{$sample[$j]}{$pos}{PROVP}=$temp[36];
		$mut{$sample[$j]}{$pos}{GERP}=$temp[49];
		$mut{$sample[$j]}{$pos}{DEPTH}=$geno[2];
		$mut{$sample[$j]}{$pos}{COSMIC}=$temp[55];
		if ($geno[0]=~/0\/1/){
		    my @count =split(/,/, $geno[1]);
		    my $af = $count[1]/$geno[2];
		    $mut{$sample[$j]}{$pos}{AF}=$af;
		    $mut{$sample[$j]}{$pos}{GENO}="1";
		}
		else{
                    $mut{$sample[$j]}{$pos}{AF}="1";
                    $mut{$sample[$j]}{$pos}{GENO}="2";
		}
	    }
	}
    }
}

#print Dumper \%mut;


my $out_file = "/data/strickt2/3102/3102.by.sample.indel.snpeff.annovar.102115.haplo.genotypes.hg19_multianno.txt";
open(OUT, ">$out_file")||die "Can't open $out_file!";

print OUT "sample\tpos\tgene\tSNV\tAA_change\tEVS6500\t1KG\tdb138\tExAC\tSIFT\tSIFT_pred\tPPHEN_HVAR\tPPHEN_HVAR_pred\tMTaster\tMTaster_pred\tMAssessor\tMAssessor_pred\tPROV\tPROV_pred\tGERP\tCOSMIC\tDEPTH\tAF\tGENO\n";
foreach my $foo (keys %mut){
    foreach my $bar (keys %{$mut{$foo}}){
	print OUT "$foo\t$bar\t$mut{$foo}{$bar}{GENE}\t$mut{$foo}{$bar}{SNV}\t$mut{$foo}{$bar}{AA}\t$mut{$foo}{$bar}{EVS6500}\t$mut{$foo}{$bar}{KG1}\t$mut{$foo}{$bar}{DB138}\t$mut{$foo}{$bar}{EXAC}\t$mut{$foo}{$bar}{SIFT}\t$mut{$foo}{$bar}{SIFTP}\t$mut{$foo}{$bar}{HVAR}\t$mut{$foo}{$bar}{HVARP}\t$mut{$foo}{$bar}{TASTE}\t$mut{$foo}{$bar}{TASTEP}\t$mut{$foo}{$bar}{MUTA}\t$mut{$foo}{$bar}{MUTAP}\t$mut{$foo}{$bar}{PROV}\t$mut{$foo}{$bar}{PROVP}\t$mut{$foo}{$bar}{GERP}\t$mut{$foo}{$bar}{COSMIC}\t$mut{$foo}{$bar}{DEPTH}\t$mut{$foo}{$bar}{AF}\t$mut{$foo}{$bar}{GENO}\n";
    }
}
