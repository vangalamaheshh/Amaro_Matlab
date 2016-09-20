# usage: perl targcovpct.pl targ.txt cov.cbb threshold outfile.txt
#
# targ.txt is a file with one line per bp in that chromosome
#      0 if not targeted
#      1 if targeted
#
# cov.cbb is a file with one line per bp in the same chromosome
#      number of reads overlapping that position
#
# threshold is an integer, e.g. 20
#
# scans through the files in synchrony and counts:
#     targ = the number of positions targeted
#     targcov = the number of targeted positions with coverage of >= threshold
#     targcovpct = 100 * targ / targcov
#
# writes to outfile.txt:
# <targ>    <targcov>     <targcovpct>%
#
# Mike Lawrence 2009-08-13

$targfile = $ARGV[0];
$covfile = $ARGV[1];
$threshold = $ARGV[2];
$outfile = $ARGV[3];
$tmpfile = $outfile . "_partial";

die "targfile $targfile not found" if (!-e $targfile);
die "covfile $covfile not found" if (!-e $covfile);
die "threshold $threshold not valid" if ($threshold<0 || $threshold>1000000000);

open(TARG,"<$targfile");
open(COV,"<$covfile");
open(OUT,">$tmpfile");

$idx = 0;
$targ = 0;
$targcov = 0;
while(<TARG>) {
    chomp;
    if (eof(COV)) { $c=0; } else { $c = <COV>; chomp $c; }
    if ($_) {
	$targ++;
	$targcov++ if ($c >= $threshold);
    }
#    print "$idx $targ $targcov\n" if (!(++$idx % 1000000));
}
close TARG;
close COV;

$targcovpct = 100 * $targcov / $targ;

printf(OUT "%d\t%d\t%.2f%%\n",$targ,$targcov,$targcovpct);
close(OUT);
rename $tmpfile, $outfile;

