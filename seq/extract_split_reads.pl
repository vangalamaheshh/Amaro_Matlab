#!/util/bin/perl -w
use strict;
my $start=$ARGV[0];
my $finish=$ARGV[1];
my $bam_file=$ARGV[2];
my $helper_fn="splitreads.helper";
my $fastq_fn="splitreads.fastq";
my $line;
my @params;

chomp($bam_file);
for (my $i=$start;$i<=$finish;$i++) {
  if (-e "$i/$helper_fn") {
    chdir $i;
    open(HELPER_FILE, "< $helper_fn");
    $line = <HELPER_FILE>;
    close(HELPER_FILE);
    @params=split(/\t/,$line);
    system("java -classpath /xchip/cga1/ydrier/assemble/GrabSplitReads:/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq:/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar GrabSplitReads $bam_file none ".$params[1]." ".(int($params[2])-int($params[9]))." ".int(int($params[2])+int($params[10]))." $fastq_fn.1\n");
    system("java -classpath /xchip/cga1/ydrier/assemble/GrabSplitReads:/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq:/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar GrabSplitReads $bam_file none ".$params[4]." ".(int($params[5])-int($params[9]))." ".int(int($params[5])+int($params[10]))." $fastq_fn.2\n");    
    system("cat $fastq_fn.1 $fastq_fn.2 > $fastq_fn");    
    chdir "..";
  }
}
