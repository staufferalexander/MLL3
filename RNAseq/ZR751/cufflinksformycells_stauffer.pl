#!/usr/bin/perl                                                                                                                                                                       
use warnings;
use strict;
use Data::Dumper;

unless ($#ARGV == 0) {print "Usage: tophat.fusion.pl <path> \n";exit;}

my $data_dir = $ARGV[0];### /data/strickt2/MLL3project/MCF12A/  /data/strickt2/3040/2018redoforstrandedness/
#my $scripts_dir = "/home/stauffkm/scripts/jobarray/";
my $main_out = "/scratch/stauffkm/MLL32018/cufflinks.refseq/";
my $GTF = "/data/strickt2/SEQreference/gtf/refgene.fixed.gtf";
my $ref = "/data/strickt2/SEQreference/hg19/hg19.complete.fa";
my $mask = "/data/strickt2/SEQreference/bundle/sorted.rrna.trna.gtf";
my $type  = "fr-firststrand";

#my @files = `find $ARGV[0] -name 'accepted_hits.bam'`; ###will this still be accepted_hits.bam, or will it be Aligned.sortedByCoord.out.bam?

#print Dumper \@files;
#my %fin;
#foreach my $i (0..$#files){
#    chomp $files[$i];
chomp $ARGV[0];
my @name = split('\/', $ARGV[0]);
my $sample = $name[5];
#/data/strickt2/MLL32016/ZR751/213/accepted_hits.bam #/data/strickt2/3040/2018redoforstrandedness/3040-TPS-216/Aligned.sortedByCoord.out.bam
#foreach my $k (0..$#name){
#$fin{$name[5]} = $ARGV[0];
#    }
#}
#print Dumper \%fin;                                                                                                                                                                        


#foreach my $foo (keys %fin){
my $out_dir = $main_out.$sample;
#    my $output_script = $scripts_dir.$foo.".cufflinks.sh";
#    open (OUT, ">$output_script")||die "Can't open $output_script!";
my $email = "kimberly.m.stauffer\@vanderbilt.edu";
#    print OUT "\#\!\/bin\/bash\n\#SBATCH --mail-user=$email\n\#SBATCH --mail-type=ALL\n\#SBATCH --nodes=1\n\#SBATCH --ntasks=8\n\#SBATCH --mem=50000M\n\#SBATCH --time=48:00:00\n\#SBATCH -o /home/stauffkm/log/$foo.cufflinks.log\n";
#`\#\!\/bin\/bash\n\#SBATCH --mail-user=$email\n\#SBATCH --mail-type=ALL\n\#SBATCH --nodes=1\n\#SBATCH --ntasks=8\n\#SBATCH --mem=50000M\n\#SBATCH --time=48:00:00\n\#SBATCH -o /home/stauffkm/log/$foo.cufflinks.log\n`;
`cufflinks -o $out_dir -G $GTF --num-threads 8 -m $mask -b $ref -u -N -q --library-type $type $ARGV[0]\n`;
#	print "cufflinks -o $out_dir -G $GTF --num-threads 8 -m $mask -b $ref -u -N  -q --library-type $type $fin{$foo}\n`;                                                                    


#	    `sbatch $output_script`;
#}
