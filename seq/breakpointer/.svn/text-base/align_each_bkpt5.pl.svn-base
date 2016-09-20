#!/util/bin/perl -w
# 
# Yotam Drier, yotamd@gmail.com
# 
use strict;
my $start=$ARGV[0];
my $finish=$ARGV[1];
my $fastaname=$ARGV[2];
my $suffix=$ARGV[3];
my $splitreads=$ARGV[4];
my $libdir=$ARGV[5];
my $bwacmd=$libdir."bwa";
my $samcmd=$libdir."samtools";

system("chmod +x $bwacmd");
system("chmod +x $samcmd");

for (my $i=$start;$i<=$finish;$i++) {
  if (-e "$i/$fastaname") {
    chdir $i;
    system("$bwacmd index $fastaname");
    system("$bwacmd aln -k 3 -n 15 $fastaname $splitreads > bwa.align.$suffix.sai");
    system("$bwacmd samse $fastaname bwa.align.$suffix.sai $splitreads > bwa.align.$suffix.sam");
    if (-s "bwa.align.$suffix.sam") {
      system("$samcmd faidx $fastaname");
      system("$samcmd import $fastaname.fai bwa.align.$suffix.sam bwa.align.$suffix.bam");
      system("$samcmd sort bwa.align.$suffix.bam bwa.align.$suffix.sorted");
      system("$samcmd view bwa.align.$suffix.sorted.bam | awk \'{if (\$6 != \"*\") print}\' > bwa.matched.$suffix.sam");
    }
    chdir "..";
  }
}
