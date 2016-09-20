#!/util/bin/perl -w
use strict;
my $start=$ARGV[0];
my $finish=$ARGV[1];
my $bamlist="";
my $line;

for (my $i=$start;$i<=$finish;$i++) {
  if (-e "$i/fused_seq.fasta") {
    chdir $i;
    system("bwa index fused_seq.fasta");
    system("bwa aln fused_seq.fasta splitreads.fastq > bwa.align.sai");
    system("bwa samse fused_seq.fasta bwa.align.sai ./splitreads.fastq > bwa.align.sam");
    if (-s "bwa.align.sam") {
      system("samtools faidx fused_seq.fasta");
      system("samtools import fused.fasta.fai bwa.align.sam bwa.align.bam");
      system("samtools sort bwa.align.bam bwa.align.sorted");
      system('samtools view bwa.align.sorted.bam | awk \'{if ($6 != "*") print}\' > bwa.matched.sam');
    }
    chdir "..";
  }
}
