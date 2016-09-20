$|=1;

# COMMANDLINE

if ($#ARGV+1!=7) {
    print "Usage: extract_from_wig.pl <targfile> <covwig> <categdir> <mincateg> <maxcateg> <outstem> <chr>\n";
    print "\n";
    print "  <targfile> should have the following structure (tab- or space-delimited)\n";
    print "        <targname/genename>   <chr>   <start>   <end>\n";
    print "\n";
    print "  <categdir> should be 'none' if no category breakdown is desired.\n";
    exit;
}

$TARGS = $ARGV[0];
$COVWIG = $ARGV[1];
$CATEGDIR = $ARGV[2];
$MINCATEG = $ARGV[3];
$MAXCATEG = $ARGV[4];
$OUTSTEM = $ARGV[5];
$CHR = $ARGV[6];

# allocate memory

$MAXLEN = 250000000;
@C = ((0) x $MAXLEN+1);
$C[$MAXLEN-1] = 0;

# read categories file

die "max <maxcateg> is 32767" if $MAXCATEG>32767;
if ($CATEGDIR ne "none") {
    $CATEGFILE = $CATEGDIR . "/chr" . $CHR . ".wig";
    print "Reading categories file $CATEGFILE\n";
    die "$CATEGFILE not found!\n" if !-e $CATEGFILE;
    open (CATEGFILE,"<$CATEGFILE");
    $_ = <CATEGFILE>;
    chomp;
    die "Category wiggle file should start with 'track' line" if !(m/^track/);
    $pos = 0;
    $chr = 0;
    $step = 0;
    while(<CATEGFILE>) {
	if (m/^\d/) {
	    die "wiggle file needs fixedStep line" if $pos==0 || $step==0 || $chr==0;
	    if ($chr==$CHR) {
		die "category out of range: $_" if $_<$MINCATEG || $_>$MAXCATEG;
		chomp;
		$C[$pos] = $_+0;
	    }
	    $pos += $step;
	    print "pos = $pos\n" if !($pos%1000000);
	} else {
	    die "variableStep not supported" if m/^variableStep/;
	    if (m/^fixedStep/) {
		$result = (($chr, $pos, $step) = m/^fixedStep chrom=chr(\S*) start=(\d*) step=(\d*)/);
		die "invalid fixedStep line in category wiggle file: $_" if $result!=3;
		$chr = 23 if $chr eq "X";
		$chr = 24 if $chr eq "Y";
		die "invalid chr in category wiggle file: $_" if $chr<1 || $chr>24;
		die "invalid step in category wiggle file: $_" if $step<1;
		die "invalid pos in category wiggle file: $_" if $pos<1;
	    }
	}
    }
    close CATEGFILE;
}		

# read coverage wiggle file

print "Reading coverage file $COVWIG\n";
die "$COVWIG not found!\n" if !-e $COVWIG;
open (COVWIG,"<$COVWIG");
$_ = <COVWIG>;
chomp;
die "Category wiggle file should start with 'track' line" if !m/^track/;
$pos = 0;
$chr = 0;
$step = 0;
while(<COVWIG>) {
    if (m/^\d/) {
	die "wiggle file needs fixedStep line" if $pos==0 || $step==0 || $chr==0;
        if ($chr==$CHR) {
	    chomp;
	    if ($_ eq "1") {
		$C[$pos] = $C[$pos] | 32768;
	    } elsif ($_ ne "0") {
		die "invalid line in coverage wiggle file: $_";
	    }
	}
        $pos += $step;
        print "pos = $pos\n" if !($pos%1000000);
    } else {
	die "variableStep not supported" if m/^variableStep/;
	if (m/^fixedStep/) {
            $result = (($chr, $pos, $step) = m/^fixedStep chrom=chr(\S*) start=(\d*) step=(\d*)/);
            die "invalid fixedStep line in coverage wiggle file: $_" if $result!=3;
            $chr = 23 if $chr eq "X";
            $chr = 24 if $chr eq "Y";
            die "invalid chr in coverage wiggle file: $_" if $chr<1 || $chr>24;
            die "invalid step in coverage wiggle file: $_" if $step<1;
            die "invalid pos in coverage wiggle file: $_" if $pos<1;
	}
    }
}
close COVWIG;

# read targets file and output totals

print "Reading targets file $TARGS and outputting totals\n";
die "$TARGS not found!" if !-e $TARGS;
open(TARGS,"<$TARGS");

$OUT = $OUTSTEM . "_chr" . $CHR;
$OUT_TEMP = $OUT . "_partial";
open(OUT,">$OUT_TEMP");

while (<TARGS>) {
    $result = (($gene,$chr,$start,$end) = m/^([^\s]*)\s*(\d*)\s*(\d*)\s*(\d*)\s*/);
    die "Invalid target: $_" if $result<4;
    $chr = 23 if $chr eq "X";
    $chr = 24 if $chr eq "Y";
    next if $chr!=$CHR;   # report only targets on the desired chromosome
    print OUT "$gene\t$chr\t$start\t$end";
    if ($CATEGDIR ne "none") {
	@T = ((0) x $MAXCATEG+1);
	for ($c=$MINCATEG;$c<=$MAXCATEG;$c++) {
          $T[$c]=0;
        }
	for ($i=$start;$i<=$end;$i++) {
	    $T[$C[$i]&32767]++ if ($C[$i]>32767);
	}
	for ($c=$MINCATEG;$c<=$MAXCATEG;$c++) {
	    print OUT "\t$T[$c]";
	}
    } else {
	$T = 0;
	for ($i=$start;$i<=$end;$i++) {
            $T++ if ($C[$i]>32767);
        }
	print OUT "\t$T";
    }
    print OUT "\n";

}
close TARGS;
close OUT;
rename $OUT_TEMP, $OUT;

