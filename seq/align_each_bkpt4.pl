#!/util/bin/perl -w
use strict;
my $start=$ARGV[0];
my $finish=$ARGV[1];
my $fastaname=$ARGV[2];
my $suffix=$ARGV[3];
my $splitreads=$ARGV[4];
my $bwacmd="/broad/software/free/Linux/suse_10_x86_64/pkgs/bwa_0.5.5/bwa";
my $samcmd="/broad/tools/Linux/x86_64/pkgs/samtools_0.1.7a/bin/samtools";

for (my $i=$start;$i<=$finish;$i++) {
  if (-e "$i/$fastaname") {
    chdir $i;
    # print("$bwacmd index $fastaname");
    system("$bwacmd index $fastaname");
    # print("$bwacmd aln -k 3 -n 15 $fastaname $splitreads > bwa.align.$suffix.sai");
    system("$bwacmd aln -k 3 -n 15 $fastaname $splitreads > bwa.align.$suffix.sai");
    # print("$bwacmd samse $fastaname bwa.align.$suffix.sai $splitreads > bwa.align.$suffix.sam\n");
    system("$bwacmd samse $fastaname bwa.align.$suffix.sai $splitreads > bwa.align.$suffix.sam");
    if (-s "bwa.align.$suffix.sam") {
      # print("$samcmd faidx $fastaname");
      system("$samcmd faidx $fastaname");
      # print("$samcmd import $fastaname.fai bwa.align.$suffix.sam bwa.align.$suffix.bam\n");
      system("$samcmd import $fastaname.fai bwa.align.$suffix.sam bwa.align.$suffix.bam");
      # print("$samcmd sort bwa.align.$suffix.bam bwa.align.$suffix.sorted");
      system("$samcmd sort bwa.align.$suffix.bam bwa.align.$suffix.sorted");
      # print("$samcmd view bwa.align.$suffix.sorted.bam | awk \'{if (\$6 != \"*\") print}\' > bwa.matched.$suffix.sam");
      system("$samcmd view bwa.align.$suffix.sorted.bam | awk \'{if (\$6 != \"*\") print}\' > bwa.matched.$suffix.sam");
    }
    chdir "..";
  }
}
