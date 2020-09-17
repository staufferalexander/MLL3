#!/usr/bin/env perl
use warnings;
use Data::Dumper;

unless ($#ARGV == 0) {print "Usage: tophat.fusion.pl <path> \n";exit;}

my $data_dir = $ARGV[0];### /data/strickt2/3040/                                                                                                                  
my $scripts_dir = "/home/stauffkm/scripts/shell/";
my $main_out = "/scratch/stauffkm/MLL32018/";
my $tmp_main = "/scratch/stauffkm/tmp/";
my $GTF = "/data/strickt2/SEQreference/gtf/refgene.fixed.gtf";  ### is this the gtf we want to use?
my $ref = "/data/strickt2/SEQreference/hg19/hg19.complete.fa";     ### is this the reference genome we want to use?                                  
my $tindex = "/data/strickt2/SEQreference/indexes/transcriptome/refseq/refgene.fixed";
my $bwt2_index = "/data/strickt2/SEQreference/indexes/bowtie2/hg19/hg19.complete"; ###was not the file that was originally used?
my $description = "fr-firststrand";
my $center = "VANTAGE";
my @files = `find $ARGV[0] -maxdepth 1 -name '*gz' `;

my $R2;
my $sample;
my $ID;
chomp $ARGV[0];
my @temp = split('\/', $ARGV[0]);
my @name2 = split(/_/, $ARGV[0]);
my @name3 = split(/_/, $temp[$#temp]);
#my $out="/".$temp[1]."/".$temp[2]."/".$temp[3]."/".$temp[4]; #/data/strickt2/3040/redo/                                                                                           
my $out="/scratch/stauffkm/MLL32018/".$temp[4];
my $R1 = $ARGV[0]; #/data/strickt2/3040/rnaseq_H5MWMBBXX/3040-TPS-212_1_sequence.txt.gz                                                                                            
if ($name2[1] == '1'){
    $R2 = $name2[0]."_2_sequence.txt.gz"; #/data/strickt2/3040/redo/3040-TPS-205                                                                                                   
    $sample = $name2[0];
    my @IDarray =  split(/-/, $sample);
    $ID = $IDarray[2];
}
else { #/data/strickt2/3040/rnaseq _ H5MWMBBXX/3040-TPS-212 _                                                                                                                      
    $R2 = $name2[0]."_".$name2[1]."_2_sequence.txt.gz";
    $sample = $name3[0];
    my @IDarray = split(/-/, $sample);
    $ID = $IDarray[2];
}
my $out_dir = $out."/".$name3[0]."/";
my $tmp_dir = $tmp_main."/".$name3[0];
`mkdir $tmp_dir`;
my $plat_unit = "Illumina_HiSeq";
my $date = "07/31/2018";
my $platform = "Illumina"; 
#/data/strickt2/3040/rnaseq_H5FW5BBXX/SampleSheet_2_sequence.txt.gz ???        
#`mkdir $tmp_dir`;
`tophat -o $out_dir -G $GTF --transcriptome-index $tindex  --num-threads 8 --keep-fasta-order --library-type fr-firststrand --tmp-dir $tmp_dir --rg-id $ID --rg-sample $sample --rg-description $description --rg-platform-unit $plat_unit --rg-center $center --rg-date $date --rg-platform $platform $bwt2_index  $R1 $R2`
#print "tophat -o $out_dir -G $GTF --transcriptome-index $tindex  --num-threads 8 --keep-fasta-order --library-type fr-firststrand --tmp-dir $tmp_dir --rg-id $ID --rg-sample $sample --rg-description $description --rg-platform-unit $plat_unit --rg-center $center --rg-date $date --rg-platform $platform $bwt2_index  $R1 $R2\n"
#my %fin;
#foreach my $i (0..$#files){
#    chomp $files[$i];
#    my @name = split('\/', $files[$i]);
#    my $dir = "/". $name[1] ."/". $name[2] ."/". $name[3] ."/". $name[4]."/";                                                                                       
#   # my $dir = "/data/strickt2/MLL32016/";#
### /data/strickt2/3040/redo/
#    my @sample1 = split('_', $name[$#name]);
#    my $lane = $sample1[1];
#    my @samplenum = split('-', $sample1[0]);
#    my $samplenumber = $samplenum[2];
#    my $laneorigin = $name[4];
#    my @sample =  split('-', $name[$#name]);
#    if ($sample[2]=~ /_1_sequence\.txt\.gz/){
#        $fin{$samplenumber}{READ1} = $files[$i];
#    }
#    if ($sample[2]=~ /_2_sequence\.txt\.gz/){
#        $fin{$samplenumber}{READ2} = $files[$i];
#    }
#  #  $fin{$samplenumber}{DIR} = "/data/strickt2/MLL3project/".$name[4]."/".$name[5];
#    $fin{$samplenumber}{DIR} = "/scratch/stauffkm/MLL32018/";
#    $fin{$samplenumber}{ID} = $sample1[0];
#    $fin{$samplenumber}{PLATFORM} = "Illumina";
#    $fin{$samplenumber}{PU} = $name[$#name]; ### what should I use here?
#    $fin{$samplenumber}{library} = "3040"; ###/data/strickt2/3040/rnaseq_H5MWMBBXX/3040-TPS-217_2_sequence.txt.gz
#    $fin{$samplenumber}{DATE} = "07/31/2018";
#    $fin{$samplenumber}{SAMPLE} = "3040_".$sample1[0]; ###
#    $fin{$samplenumber}{CENTER} = "VANTAGE";
#    $fin{$samplenumber}{LANE} = $lane;
#    $fin{$samplenumber}{DESCRIPTION} = "MLL32018_RNASeq";
#}
#
#print Dumper \%fin;
#
#foreach my $foo (keys %fin){
#    my $out_dir = $main_out.$foo;
#    mkdir $out_dir;
##    print "$out_dir\n";
#    my $tmp_dir = $tmp_main.$foo;
#    mkdir $tmp_dir;
##    print "$tmp_dir\n";
#    my $output_script = $scripts_dir.$foo."tophat.fusion.sh";
##    print "$output_script\n";
#    open (OUT, ">$output_script")||die "Can't open $output_script!";
#    my $email = "kimberly.m.stauffer\@vanderbilt.edu";
#    print OUT "\#\!\/bin\/bash\n\#SBATCH --mail-user=$email\n\#SBATCH --mail-type=ALL\n\#SBATCH --nodes=1\n\#SBATCH --ntasks=8\n\#SBATCH --mem=50000M\n\#SBATCH --time=48:00:00\n\##SBATCH -o /home/stauffkm/log/$foo.tophat.log\n";
#  #  print "\#\!\/bin\/bash\n\#PBS -M $email\n\#PBS -m bae\n\#PBS -l nodes=1:ppn=8\n\#PBS -l mem=30000mb\n\#PBS -l walltime=24:00:00\n\#PBS -o /home/stauffkm/log/$foo.tophat.log\n\\#PBS -j oe\n";                                                                                                                                                                     
#    print OUT  "tophat -o $out_dir -G $GTF --transcriptome-index $tindex  --num-threads 8 --keep-fasta-order --library-type fr-firststrand --tmp-dir $tmp_dir --rg-id $fin{$foo}{ID} --rg-sample $fin{$foo}{SAMPLE} --rg-description $fin{$foo}{DESCRIPTION} --rg-platform-unit $fin{$foo}{PU} --rg-center $fin{$foo}{CENTER} --rg-date $fin{$foo}{DATE} --rg-platform $fin{$foo}{PLATFORM} $bwt2_index  $fin{$foo}{READ1} $fin{$foo}{READ2}\n";
#    print "tophat -o $out_dir -G $GTF --transcriptome-index $tindex  --num-threads 8 --keep-fasta-order --library-type fr-firststrand --tmp-dir $tmp_dir --rg-id $fin{$foo}{ID} --rg-sample $fin{$foo}{SAMPLE} --rg-description $fin{$foo}{DESCRIPTION} --rg-platform-unit $fin{$foo}{PU} --rg-center $fin{$foo}{CENTER} --rg-date $fin{$foo}{DATE} --rg-platform $fin{$foo}{PLATFORM} $bwt2_index  $fin{$foo}{READ1} $fin{$foo}{READ2}\n\n\n";
##    print OUT "perl /home/stauffkm/scripts/remove.pl $out_dir\n";
#   # print OUT "perl /home/stauffkm/scripts/toremovelogs.pl $out_dir\n";

#    close OUT;
#    `sbatch $output_script`;
#}

