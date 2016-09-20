#!/util/bin/perl -w
use strict;
my $start=$ARGV[0];
my $finish=$ARGV[1];
my $bam_file=$ARGV[2];
my $helper_fn="splitreads.helper2";
my $fastq_fn="splitreads2.fastq";
my ($line, $take_before_1, $take_before_2, $take_after_1, $take_after_2);
my @params;

chomp($bam_file);
for (my $i=$start;$i<=$finish;$i++) {
  if (-e "$i/$helper_fn") {
    chdir $i;
    open(HELPER_FILE, "< $helper_fn");
    $line = <HELPER_FILE>;
    close(HELPER_FILE);
    @params=split(/\t/,$line);    
    if ($params[3]==0) {
    	$take_before_1=int($params[9])+int($params[10]);
    	$take_after_1=int($params[9]);
    } else {
    	$take_before_1=int($params[9]);
    	$take_after_1=int($params[9])+int($params[10]);
    }
    if ($params[6]==0) {
    	$take_before_2=int($params[9])+int($params[10]);
       	$take_after_2=int($params[9]);
    } else {
        $take_before_2=int($params[9]);
        $take_after_2=int($params[9])+int($params[10]);
    }
    system("java -classpath /xchip/cga1/ydrier/assemble/GrabSplitReads:/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq:/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar GrabSplitReads $bam_file none ".$params[1]." ".(int($params[2])-$take_before_1)." ".int(int($params[2])+$take_after_1)." $fastq_fn.1\n");
    system("java -classpath /xchip/cga1/ydrier/assemble/GrabSplitReads:/home/radon00/lawrence/CancerGenomeAnalysis/trunk/matlab/seq:/xchip/tcga/gbm/analysis/lawrence/samtools/sam-1.07.jar GrabSplitReads $bam_file none ".$params[4]." ".(int($params[5])-$take_before_2)." ".int(int($params[5])+$take_after_2)." $fastq_fn.2\n");    
    system("cat $fastq_fn.1 $fastq_fn.2 > $fastq_fn");    
    chdir "..";
  }
}
