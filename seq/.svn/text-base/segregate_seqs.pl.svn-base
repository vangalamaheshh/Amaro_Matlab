# usage: perl segregate_seqs infile numsoutfile seqsoutfile
#
# infile should be "all.weird.pairs", with 13 numeric cols followed by 2 string cols
# numsoutfile will contain just the numeric cols
# seqsoutfile will contain just the string cols
#
# Mike Lawrence 2009-07-31

$IFILE = $ARGV[0];
$NFILE = $ARGV[1];
$SFILE = $ARGV[2];
$NFILE_TEMP = $NFILE . "_partial";
$SFILE_TEMP = $SFILE . "_partial";

if (!(-e $IFILE)) { print "$INFILE not found!\n"; exit; }
open(IN,"<$IFILE");
open(NOUT,">$NFILE_TEMP");
open(SOUT,">$SFILE_TEMP");

$line = 0 ;
$pattern = "(" . "[^\t]*\t"x12 . "[^\t]*)\t([^\t]*\t[^\t]*)";
while(<IN>) {
    $line++;
    chomp;
    $result = (($n,$s) = m/^$pattern$/);
    if (!$result) { print "$_\nmatch failed at line $line\n"; exit; }
    print NOUT "$n\n";
    print SOUT "$s\n";
}

close IN;
close NOUT;
close SOUT;
rename $NFILE_TEMP, $NFILE;
rename $SFILE_TEMP, $SFILE;
