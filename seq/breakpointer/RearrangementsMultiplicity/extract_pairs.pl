#!/usr/bin/perl -w
use strict;

#my $bam="/xchip/cga1/firehose_output/trunk/Individual/".$sample."/wgs/bam/raw/tumor.bam";
die "Usage: extract_pairs <bam file name> <rearrangements file> <blacklist file> <window size> [libdir]\nExtract pairs around rearrngements (+-window size).\nRearrangements file should contain header line and columns <chr1> <pos1> <chr2> <pos2> in that order.\nOutputs under current directory.\nNot enough parameters" if ($#ARGV+1 < 4);
my $bam=$ARGV[0];
my $filename=$ARGV[1];
my $blacklist=$ARGV[2];
my $expand_pairs_extraction = $ARGV[3];
my $libdir = "";
if ($#ARGV+1 >= 5) {
    $libdir = $ARGV[4];
}

#my $filename="/xchip/tcga_scratch/lawrence/pr/".$sample."/wgs/dRanger_results.txt";
#my $cmd="java -classpath /home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq:/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar GrabPairs /xchip/tcga_scratch/lawrence/pr/".$sample."/wgs/tumor.bam $blacklist";
my $cmd="java -classpath /xchip/cga1/ydrier/CancerGenomeAnalysis/trunk/analysis_pipeline/tools/dist/GrabSplitReads.jar:/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar GrabSplitReads $bam $blacklist";

open(INFILE, "< $filename");
my $i=0;
my $header=<INFILE>;
foreach my $thisline (<INFILE>) {
  my @this_line = split("\t", $thisline);
  $i++;
  mkdir("$i",0755);
  my $chr1=$this_line[0];
  my $chr2=$this_line[2];
  my $beg1=$this_line[1]-$expand_pairs_extraction;
  my $end1=$this_line[1]+$expand_pairs_extraction;
  my $beg2=$this_line[3]-$expand_pairs_extraction;
  my $end2=$this_line[3]+$expand_pairs_extraction;
  print("$thisline\n");
  system("$cmd $chr1 $beg1 $end1 $i/first.txt");
 #   print("$cmd $chr1 $beg1 $end1 $i/first.txt\n");
  system("$cmd $chr2 $beg2 $end2 $i/second.txt");
 #   print("$cmd $chr2 $beg2 $end2 $i/second.txt\n");
  system("cat $i/first.txt $i/second.txt | sort | perl /home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq/joindump.pl 1 > $i/bothends.txt");
}
close(INFILE);
print "$i\n"; 
  
