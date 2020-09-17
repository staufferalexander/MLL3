#!/usr/bin/perl
use warnings;
use Data::Dumper;

unless ($#ARGV == 0) {print "Usage: bwa.pl < path >\n";exit;} ### /scratch/stauffkm/chipseq/
#/data/strickt2/MLL32501/
my $data_dir = $ARGV[0];
my $out_dir = $data_dir."alignments/";  ### data_dir is /data/strickt2/chipseq/chipseq_H5CC5BBXX/
my $scripts_dir = "/home/stauffkm/scripts/shell/";

my $ref_fa = "/data/strickt2/SEQreference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa";
my $ref_fai = "/data/strickt2/SEQreference/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa.fai";

my @files = `find $ARGV[0] -name '*10k.fa'`;
#$my @files = `find $ARGV[0] -name '*.gz'`;
#my @files = ("/data/strickt2/chipseq/HFNH7BBXX/3040-TPS-191_1_sequence.txt.gz", "/data/strickt2/chipseq/HFNH7BBXX/3040-TPS-137_1_sequence.txt.gz", "/data/strickt2/chipseq/HFNH7BBXX/3040-TPS-134_1_sequence.txt.gz", "/data/strickt2/chipseq/HFNH7BBXX/3040-TPS-133_1_sequence.txt.gz", "/data/strickt2/chipseq/HFNH7BBXX/3040-TPS-129_1_sequence.txt.gz");

my %fin;
foreach my $i (0..$#files){
    chomp $files[$i];
    my @name = split('\/', $files[$i]);
    my @sample = split('\_', $name[$#name]);
    $fin{$sample[0]} = $files[$i];
}

print Dumper \%fin;
#`mkdir $out_dir`;
foreach my $foo (keys %fin){
    my $output_script = $scripts_dir.$foo.".bwa.sh";
    open (OUT, ">$output_script")||die "Can't open $output_script!";
    my $email = "kimberly.m.stauffer\@vanderbilt.edu";
    my $SE_out = $data_dir."alignments/".$foo.".SE.sam";
    my $SE_out_bam = $data_dir."alignments/".$foo.".SE.bam";
#    print "$SE_out\t$SE_out_bam\n";
    print OUT "\#\!\/bin\/bash\n\#SBATCH --mail-user=$email\n\#SBATCH --mail-type=ALL\n\#SBATCH --nodes=1\n\#SBATCH --ntasks=8\n\#SBATCH --mem=120000mb\n\#SBATCH --time=6:00:00\n\#SBATCH -o /home/stauffkm/log/$foo.bam.log\n";
    print OUT  "bwa mem -t 8 -M $ref_fa $fin{$foo} > $SE_out\n";
#    print  "bwa mem -t 8 -M $ref_fa $fin{$foo} > $SE_out\n";
    print OUT "samtools view -bSt $ref_fai $SE_out -o $SE_out_bam\n";
    print OUT "rm $SE_out\n";
#    `sbatch $output_script`;
}
