# COMMANDLINE

if ($#ARGV+1!=7) {
    print "Usage: get_TN_coverage_from_cbb.pl <targfile> <tumdir_cbb> <normdir_cbb> <categdir> <numcategs> <outstem> <chr>\n";
    print "\n";
    print "  <targfile> should have the following structure (tab- or space-delimited)\n";
    print "        <targname/genename>   <chr>   <start>   <end>\n";
    print "\n";
    print "  <categdir> should be 'none' if no category breakdown is desired.\n";
    exit;
}

$TARGS = $ARGV[0];
$TUMDIR = $ARGV[1];
$NORMDIR = $ARGV[2];
$CATEGDIR = $ARGV[3];
$NUMCATEGS = $ARGV[4];
$OUTSTEM = $ARGV[5];
$CHR = $ARGV[6];

if (!(-e $TARGS)) { print "$TARGS not found!\n"; exit; }
open(TARGS,"<$TARGS");
$i=0;
while (<TARGS>) {
    $result = (($gene[$i],$chr,$start[$i],$end[$i]) = m/^([^\s]*)\s*(\d*)\s*(\d*)\s*(\d*)\s*/);
    if (!$result) { print "Invalid target: $_\n"; exit; }
    if ($chr==$CHR) { $i++; }   # only keep targets on the desired chromosome
}
close TARGS;
$nt = $i;
if (!$nt) { print "No targets on chromosome $CHR!\n"; exit; }

# make list of positions queried
print "Analyzing targets...\n";
$last_pos_queried = 0;
$collisions = 0;
for ($i=0;$i<$nt;$i++) {
    for ($c=0;$c<$NUMCATEGS;$c++) {
      $NORMcov[$i][$c] = 0;   # initialize counts
      $TUMcov[$i][$c] = 0;
    }
    for ($j=$start[$i];$j<=$end[$i];$j++) {
	if ($pos_queried[$j]>0) { $collisions++; }   # this position was already assigned to a target
	else { $pos_queried[$j] = $i+1; }   # use one-based index, so zero can indicate "none"
	if ($j>$last_pos_queried) { $last_pos_queried = $j; }
    }
}
if ($collisions) {
  print "Warning: $collisions bases were specified in more than one target.\n".
    "For these positions, coverage will be credited to only the first target specified.\n";
}

# CATEGORIES

if ($CATEGDIR eq "none") {
    if ($NUMCATEGS!=1) {
	print "Warning: <numcategs> of $NUMCATEGS is meaningless without specifying <categdir>";
        $NUMCATEGS = 1;
    }
} else {
    if ($NUMCATEGS<2) {
	print "Warning: <numcategs> of $NUMCATEGS is meaningless when specifying <categidr>";
    }
    $CATEGFILE = $CATEGDIR . "/chr" . $CHR . ".txt";
    if (!(-e $CATEGFILE)) { print "$CATEGFILE not found!\n"; exit; }
    open (CATEGFILE,"<$CATEGFILE");
}

# SCAN COVERAGE

$TUMFILE =   $TUMDIR . "/chr" . $CHR . ".cbb";
$NORMFILE = $NORMDIR . "/chr" . $CHR . ".cbb";

if (!(-e $TUMFILE)) { print "$TUMFILE not found!\n"; exit; }
if (!(-e $NORMFILE)) { print "$NORMFILE not found!\n"; exit; }

open (TUMFILE,"<$TUMFILE");
open (NORMFILE,"<$NORMFILE");

$pos = 0;
$categ = 1;
while ($pos<=$last_pos_queried) {
    $tum = <TUMFILE>;
    $norm = <NORMFILE>;
    if ($NUMCATEGS>1) {
	$categ = <CATEGFILE>;
    }
    $pos++;
    if (!($pos%1000000)) { print "chr$CHR:$pos\n"; }
    $i = $pos_queried[$pos]-1;
    if ($i==-1) { next; }  # base not in a targeted region
    chomp $tum;
    chomp $norm;
    chomp $categ;
    $TUMcov[$i][$categ-1] += $tum;
    $NORMcov[$i][$categ-1] += $norm;
}

close TUMFILE;
close NORMFILE;
if ($CATEGDIR ne "none") { close CATEGFILE; }

# write report

$OUT = $OUTSTEM . "_chr" . $CHR;
$OUT_TEMP = $OUT . "_partial";
open(OUT,">$OUT_TEMP");
for ($i=0;$i<$nt;$i++) {
    print OUT "$gene[$i]\t$CHR\t$start[$i]\t$end[$i]";
    for ($c=0;$c<$NUMCATEGS;$c++) {
	print OUT "\t$TUMcov[$i][$c]";
    }
    for ($c=0;$c<$NUMCATEGS;$c++) {
        print OUT "\t$NORMcov[$i][$c]";
    }
    print OUT "\n";
}
close OUT;
rename $OUT_TEMP, $OUT;

