#!/util/bin/perl -w
use strict;
my $line;
my $matches;
my $start=$ARGV[0];
my $finish=$ARGV[1];
my $allfused_fn="all_fused.fasta";
my $sam_fn="all_aligned.sam";
my $matches_fn="count_matches.txt";

open(ALLFUSED_FILE, "> $allfused_fn");
open(SAMFILE, "> $sam_fn");
open(COUNTMATCHES, "> $matches_fn");
for (my $i=$start;$i<=$finish;$i++) {
  if (-e "$i/fused_seq.fasta") {
    open(INFILE, "< $i/fused_seq.fasta");
    foreach $line (<INFILE>) {
      print ALLFUSED_FILE $line;
    }
    close(INFILE);
    if (-s "$i/bwa.align.sam") {
      open(INFILE, "< $i/bwa.matched.sam");
      foreach $line (<INFILE>) {
	print SAMFILE $line;
      }
      close(INFILE);
      $matches = `cat $i/bwa.matched.sam | wc -l`;
      print COUNTMATCHES "$i\t$matches";
    }
  }
}
close(ALLFUSED_FILE);
close(SAMFILE);
close(COUNTMATCHES);
system("samtools faidx $allfused_fn");
system("samtools import ./$allfused_fn.fai ./$sam_fn ./bwa.matched.bam");
system("samtools index bwa.matched.bam");
