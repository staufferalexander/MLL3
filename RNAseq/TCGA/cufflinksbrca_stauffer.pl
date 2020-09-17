#!/usr/bin/perl                                                                                                                                                                                       
use warnings;
use strict;
use Data::Dumper;

unless ($#ARGV == 0) {print "Usage: cufflinksbrca.pl <path> \n";exit;} ### /scratch/strickt2/tcga/

my $data_dir = $ARGV[0];
my $scripts_dir = "/home/stauffkm/scripts/shell/";
my $main_out = "$ARGV[0]";
my $GTF = "/data/strickt2/SEQreference/bundle/gtf/hg19.complete.new.gtf";
my $ref = "/data/strickt2/SEQreference/indexes/bowtie2/hg19/hg19.complete.fa";
my $mask = "/data/strickt2/SEQreference/bundle/sorted.rrna.trna.gtf";
my $type  = "fr-unstranded";

my @files = `find $ARGV[0] -name 'accepted_hits.bam'`;

print Dumper \@files;

my %fin;
foreach my $i (0..$#files){
    chomp $files[$i];
    my @name = split('\/', $files[$i]);
    $fin{$name[4]} = $files[$i]; 
}

foreach my $foo (keys %fin){
    my $out_dir = $main_out.$foo;                                   
    my $output_script = $scripts_dir.$foo.".cufflinks.sh";
    open (OUT, ">$output_script")||die "Can't open $output_script!";
    my $email = "kimberly.m.stauffer\@vanderbilt.edu";
    print OUT "\#\!\/bin\/bash\n\#PBS -M $email\n\#PBS -m bae\n\#PBS -l nodes=1:ppn=8\n\#PBS -l mem=20000mb\n\#PBS -l walltime=12:00:00\n\#PBS -o /home/stauffkm/log/$foo.cufflinks.log\n\#PBS -j oe\n";
    print OUT  "cufflinks -o $out_dir -G $GTF --num-threads 8 -m $mask -b $ref -u -q --library-type $type $fin{$foo}\n";
    print "cufflinks -o $out_dir -G $GTF --num-threads 8 -m $mask -b $ref -u -q --library-type $type $fin{$foo}\n";

    #`qsub $output_script`;                                                                                                                                                                           
}
