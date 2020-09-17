#!/usr/bin/perl                                                                                                                                                                           
use warnings;
use Data::Dumper;

unless ($#ARGV == 0) {print "Usage: tophat.fusion.pl <path> \n";exit;}

my $data_dir = $ARGV[0];### /scratch/strickt2/tcga/
my $scripts_dir = "/home/stauffkm/scripts/shell/";
my $main_out = $ARGV[0];
my $tmp_main = "/scratch/strickt2/tmp/";
my $GTF = "/data/strickt2/SEQreference/gtf/refgene.fixed.gtf";
my $ref = "/data/strickt2/SEQreference/indexes/bowtie2/hg19/hg19.complete.fa";                                                                                                                                                                                          
my $tindex = "/data/strickt2/SEQreference/indexes/transcriptome/refseq/";
my $bwt2_index = "/data/strickt2/SEQreference/indexes/bowtie2/hg19/hg19.complete";
my $description = "TruSeq.RNA.unstranded.polyA"; 
my $center = "UNC.CH"; 

my @files = `find $ARGV[0] -name '*fastq.gz' `;

my %fin;
foreach my $i (0..$#files){
    chomp $files[$i];
    my $lane;
    my @name = split('\/', $files[$i]);
    #print Dumper \@name; 
    my $dir = "/". $name[1] ."/". $name[2] ."/". $name[3] ."/". $name[4];                                                                                                                                                                
    my @sample =  split('_', $name[5]);
    #print Dumper \@sample;
    my $length = scalar(@sample); 
    my $li = $length - 2;
    my $li2 = $length - 1;
    if ($sample[$li] =~ /L\d+/){  
	$lane = $sample[$li]; 
    }
    elsif ($sample[$li] =~ /\./){
	my @sam = split (/\./, $sample[$li]);
	$lane =$sam[1]; 
    }
    #print Dumper \@sample;                                                                                                                                                               
    if ($sample[$li2] =~ /1\.fastq\.gz$/){  
        $fin{$name[4]}{READ1} = $files[$i];  
    }
    elsif ($sample[$li2] =~ /2\.fastq\.gz$/){   
        $fin{$name[4]}{READ2} = $files[$i];    
    }
    $fin{$name[4]}{DIR} = '/scratch/strickt2/tcga/$name[4]';
    $fin{$name[4]}{ID} = $name[4];
    $fin{$name[4]}{PLATFORM} = "Illumina";
    $fin{$name[4]}{PU} = $name[5];       
    $fin{$name[4]}{library} = $name[4];
    $fin{$name[4]}{DATE} = $sample[1];
    $fin{$name[4]}{SAMPLE} = $name[4];
    $fin{$name[4]}{CENTER} = "UNC";
    $fin{$name[4]}{LANE} = $lane;
    $fin{$name[4]}{DESCRIPTION} = "TCGA_RNASeq";
}

print Dumper \%fin;
  
foreach my $foo (keys %fin){
    my $out_dir = $main_out.$foo;
    mkdir $out_dir;
    my $tmp_dir = $tmp_main.$foo;
    mkdir $tmp_dir;
    my $output_script = $scripts_dir.$foo."tophat.fusion.sh";
    open (OUT, ">$output_script")||die "Can't open $output_script!";
    my $email = "kimberly.m.stauffer\@vanderbilt.edu";
    print OUT "\#\!\/bin\/bash\n\#PBS -M $email\n\#PBS -m bae\n\#PBS -l nodes=1:ppn=8\n\#PBS -l mem=30000mb\n\#PBS -l walltime=48:00:00\n\#PBS -o /home/stauffkm/log/$foo.tophat.log\n\#PBS -j oe\n";
  #  print "\#\!\/bin\/bash\n\#PBS -M $email\n\#PBS -m bae\n\#PBS -l nodes=1:ppn=8\n\#PBS -l mem=30000mb\n\#PBS -l walltime=24:00:00\n\#PBS -o /home/stauffkm/log/$foo.tophat.log\n\#PBS -j oe\n";
    print  OUT  "tophat -o $out_dir -G $GTF --transcriptome-index $tindex  --num-threads 8 --keep-fasta-order --library-type fr-unstranded --tmp-dir $tmp_dir --rg-id $fin{$foo}{ID} --rg-sample $fin{$foo}{SAMPLE} --rg-description $fin{$foo}{DESCRIPTION} --rg-platform-unit $fin{$foo}{PU} --rg-center $fin{$foo}{CENTER} --rg-date $fin{$foo}{DATE} --rg-platform $fin{$foo}{PLATFORM} $bwt2_index  $fin{$foo}{READ1} $fin{$foo}{READ2}\n";   
 #   print "tophat -o $out_dir -G $GTF --transcriptome-index $tindex --num-threads 8 --keep-fasta-order --library-type fr-unstranded --tmp-dir $tmp_dir --rg-id $fin{$foo}{ID} --rg-sample $fin{$foo}{SAMPLE} --rg-description $fin{$foo}{DESCRIPTION} --rg-platform-unit $fin{$foo}{PU} --rg-center $fin{$foo}{CENTER} --rg-date $fin{$foo}{DATE} --rg-platform $fin{$foo}{PLATFORM} $bwt2_index $fin{$foo}{READ1} $fin{$foo}{READ2}\n";  

    print OUT "perl /home/stauffkm/scripts/remove.pl $out_dir\n";
    print OUT "perl /home/stauffkm/scripts/toremovelogs.pl $out_dir\n";
    
    close OUT;

    `qsub $output_script`;
}

