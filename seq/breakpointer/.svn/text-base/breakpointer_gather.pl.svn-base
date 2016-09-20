#!/util/bin/perl -w
#
# Yotam Drier, yotamd@gmail.com
#
use strict;

die "Usage: breakpointer_gather <sample.id> <libdir> <resultsdir> [<keeptempfiles>]\n" if ($#ARGV+1 < 3);
my $sample = $ARGV[0];
my $libdir = $ARGV[1];
my $resultsdir = $ARGV[2];
my $keeptempfiles = 0; 
if ($#ARGV+1 >= 4) {
    $keeptempfiles = ($ARGV[3] eq "1" || $ARGV[3] eq "true");
}

my ($start, $finish, $i, $dir, $line);
my @t;
my @dirs = `ls -dX1 $resultsdir/scatter.*/`;
my $breakpoints_in_fn="my_breakpoints_all.txt";
my $allfused_fn=$sample.".fused.fasta";
my $sam_fn=$sample.".matched.sam";
my $breakpoints_out_fn=$sample.".breakpoints.txt";
my $bam_fn=$sample.".bwa.matched.bam";
my $fusedseq_fn="fused_seq.fasta";
my $seqs_inf_fn="fused_seq.inf";

open(ALLBKPTS_FILE, "> $resultsdir/$breakpoints_out_fn");
print ALLBKPTS_FILE "bpnum\tdrnum\tchr1\tpos1_orig\tchr2\tpos2_orig\tinversion\tbpreads\tqual\tseq2first\tfs_len\tfs\tpos1\tpos2\tbwareads\texact_mh\tsoft_mh\tkept1\tdropped1\tdropped2\tkept2\n";
open(ALLFUSED_FILE, "> $resultsdir/$allfused_fn");
open(SAMFILE, "> $resultsdir/$sam_fn");

if ($sample ne "NO_DATA") {
    foreach $dir (@dirs) {
	chomp($dir);
	chdir($dir);
	@t=split(/\t/,`ls -1d */ | sort -n | head -1`);
	($start) = ($t[0] =~ m/(\d+)\/$/);
	@t=split(/\t/,`ls -1d */ | sort -n | tail -1`);
	($finish) = ($t[0] =~ m/(\d+)\/$/);
	chdir("..");
	my @matches = (0) x ($finish-$start+1);
	my @soft_mh = (-1) x ($finish-$start+1);
	my @exact_mh = (-1) x ($finish-$start+1);
	my @pos1 = (-1) x ($finish-$start+1);
	my @pos2 = (-1) x ($finish-$start+1);
	my @seq_kept1 = ("failed") x ($finish-$start+1);
	my @seq_kept2 = ("failed") x ($finish-$start+1);
	my @seq_dropped1 = ("failed") x ($finish-$start+1);
	my @seq_dropped2 = ("failed") x ($finish-$start+1);
	for ($i=$start;$i<=$finish;$i++) {
	    if (-e "$dir/$i/$fusedseq_fn") {
		open(INFILE, "< $dir/$i/$fusedseq_fn");
		foreach $line (<INFILE>) {
		    print ALLFUSED_FILE $line;
		}
		close(INFILE);
		if (-s "$dir/$i/bwa.matched.bp.sam") {
		    open(INFILE, "< $dir/$i/bwa.matched.bp.sam");
		    foreach $line (<INFILE>) {
			print SAMFILE $line;
		    }
		    close(INFILE);
		    $matches[$i-$start] = `cat $dir/$i/bwa.matched.bp.sam | wc -l`;
		    chomp($matches[$i-$start]);
		}
		open(INFILE, "< $dir/$i/$seqs_inf_fn");
		foreach $line (<INFILE>) {
		    chomp($line);
		    @t=split(/\t/,$line);
		    $pos1[$i-$start] = $t[3];
		    $pos2[$i-$start] = $t[5];
		    $seq_kept1[$i-$start] = $t[6];
		    $seq_dropped1[$i-$start] = $t[7];
		    $seq_dropped2[$i-$start] = $t[9];
		    $seq_kept2[$i-$start] = $t[10];
	            $exact_mh[$i-$start] = $t[11];
		    $soft_mh[$i-$start] = $t[12];
		}
		close(INFILE);
	    }
	    system("rm -rf $dir/$i") if $keeptempfiles==0;
	}
	open(INFILE, "< $dir/$breakpoints_in_fn");
	foreach $line (<INFILE>) {
	    chomp($line);
	    @t=split(/\t/,$line);
	    $i=$t[0]-$start;
	    print ALLBKPTS_FILE "$line\t$pos1[$i]\t$pos2[$i]\t$matches[$i]\t$exact_mh[$i]\t$soft_mh[$i]\t$seq_kept1[$i]\t$seq_dropped1[$i]\t$seq_dropped2[$i]\t$seq_kept2[$i]\n";
	}
	close(INFILE);
	system("rm $dir/$breakpoints_in_fn") if $keeptempfiles==0;
    }
}

close(ALLFUSED_FILE);
close(SAMFILE);
close(ALLBKPTS_FILE);

die "No breakpoints file generated" if (!-e "$resultsdir/$breakpoints_out_fn");
die "No fused fasta generated" if (!-e "$resultsdir/$allfused_fn");
die "No sam file generated" if (!-e "$resultsdir/$sam_fn");

if ((-z "$resultsdir/$allfused_fn") or ($sample eq "NO_DATA")) { # make dummy output files
  open(TMP, "> $resultsdir/$allfused_fn.fai"); close(TMP);
  open(TMP, "> $resultsdir/$bam_fn"); close(TMP);
  open(TMP, "> $resultsdir/$bam_fn.bai"); close(TMP);
} else {
  system("chmod +x ".$libdir."samtools");
  system($libdir."samtools faidx $resultsdir/$allfused_fn");
  system($libdir."samtools import $resultsdir/$allfused_fn.fai $resultsdir/$sam_fn $resultsdir/$bam_fn 2>&1");
  system($libdir."samtools index $resultsdir/$bam_fn");   
  die "No bam generated" if (!-e "$resultsdir/$bam_fn");
} 
