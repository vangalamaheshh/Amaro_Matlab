# THRESHOLDS

$TUM_COVERAGE_NECESSARY = 14;
$NORM_COVERAGE_NECESSARY = 8;

# COMMANDLINE

if ($#ARGV+1<7) {
    print "Usage: get_coverage_stats.pl <tumdir_cbb> <normdir_cbb> <filterdir> <outdir> <chr> <mincateg> <maxcateg>\n";
    print "  <filterdir> should be a 'categdir' like\n";
    print "              /xchip/tcga_scratch/lawrence/db/context\n";
    print "              or else 'none' for no category-based filtering\n";
    print "\n";
    print "Output:  Chr  Categ   Terr   Tseqbp   Nseqbp  Callablebp\n";
    print "   (one line per category)\n";            
    print "   Terr = territory (number of bp) on this chr for this category\n";
    print "   Tseqbp = total number of sequenced bases for this chr/categ in the tumor\n";
    print "   Nseqbp = total number of sequenced bases for this chr/categ in the normal\n";
    print "   Callablebp = number of bases sequenced >=14x in tumor and >=8x in normal\n";
    exit;
}

$TUMDIR = $ARGV[0];
$NORMDIR = $ARGV[1];
$FILTERDIR = $ARGV[2];
$OUTDIR = $ARGV[3];
$CHR = $ARGV[4];
$MINCATEG = $ARGV[5];
$MAXCATEG = $ARGV[6];

if ($CHR<1 || $CHR>24) { print "chr should be between 1 and 24\n"; exit; }

# OPEN FILES

$TUMFILE =   $TUMDIR . "/chr" . $CHR . ".cbb";
if (!(-e $TUMFILE)) {
    $TUMFILE =   $TUMDIR . "/chr" . $CHR . ".txt";
    if (!(-e $TUMFILE)) {
	print "$TUMFILE not found!\n"; exit;
    }
}

$NORMFILE =   $NORMDIR . "/chr" . $CHR . ".cbb";
if (!(-e $NORMFILE)) {
    $NORMFILE =   $NORMDIR . "/chr" . $CHR . ".txt";
    if (!(-e $NORMFILE)) {
        print "$NORMFILE not found!\n"; exit;
    }
}

open (TUMFILE,"<$TUMFILE");
open (NORMFILE,"<$NORMFILE");

if ($FILTERDIR eq "none") {
    $filtering = 0;
} else {
    $filtering = 1;
    $FILTERFILE = $FILTERDIR ."/chr" .$CHR. ".txt";
    if (!(-e $FILTERFILE)) { print "$FILTERFILE not found!\n"; exit; }
    open (FILTERFILE,"<$FILTERFILE");
}

# SCAN COVERAGE

@terr = ((0) x ($MAXCATEG+1));
@tseqbp = ((0) x ($MAXCATEG+1));
@nseqbp = ((0) x ($MAXCATEG+1));
@callable = ((0) x ($MAXCATEG+1));

$pos = 0;
$categ = 1;

while (!eof(FILTERFILE)) {
    if (!eof(TUMFILE)) { $tum = <TUMFILE>; chomp $tum; } else { $tum = 0; }
    if (!eof(NORMFILE)) { $norm = <NORMFILE>; chomp $norm; } else { $norm = 0; }
    if ($filtering) { $categ = <FILTERFILE>; chomp $categ; }

    $pos++;
    if (!($pos%1000000)) { print "chr$CHR:$pos\n"; }

    $terr[$categ]++;
    $tseqbp[$categ]+=$tum;
    $nseqbp[$categ]+=$norm;
    if ($tum>=$TUM_COVERAGE_NECESSARY && $norm>=$NORM_COVERAGE_NECESSARY) {
	$callable[$categ]++;
    }
}

close TUMFILE;
close NORMFILE;
if ($filtering) { close FILTERFILE; }

# write output

$OUTFILE = $OUTDIR . "/chr" . $CHR . ".stats";
open(OUTFILE,">$OUTFILE");
for ($i=$MINCATEG;$i<=$MAXCATEG;$i++) {
    print OUTFILE "$CHR\t$i\t$terr[$i]\t$tseqbp[$i]\t$nseqbp[$i]\t$callable[$i]\n";
}
close OUTFILE;


