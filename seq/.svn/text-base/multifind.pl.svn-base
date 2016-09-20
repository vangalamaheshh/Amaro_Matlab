$TARGS = $ARGV[0];
$IN = $ARGV[1];
$OUT = $ARGV[2];
$OUT_TEMP = $OUT . "_partial";

open(TARGS,"<$TARGS");
$i=0;
while (<TARGS>) {
    $result = (($chr[$i],$start[$i],$end[$i]) = m/^\s*(\d*)\s*(\d*)\s*(\d*)$/);
    if (!$result) { print "Invalid target: $_\n"; exit; }
    $i++;
}
close TARGS;
$nt = $i;

if (!$nt) { print "No targets!\n"; exit; }
if (!(-e $IN)) { print "Input file not found!\n"; exit; }

open(OUT,">$OUT_TEMP");
open (IN,"<$IN");
while (<IN>) {
    $result = (($chr1,$start1,$end1,$chr2,$start2,$end2) = 
       m/^\d*\t(-?\d*)\t-?.\t(-?\d*)\t(-?\d*)\t(-?\d*)\t-?.\t(-?\d*)\t(-?\d*)/);
#    print "$chr1:$start1-$end1 + $chr2:$start2-$end2\n";
    if (!$result) { print $_; next; }
#    print "CHECKING $nt entries:\n";
    for ($i=0;$i<$nt;$i++) {
#	print "    ($i) $chr[$i]:$start[$i]-$end[$i]\n";
        if (($chr1==$chr[$i] && $start1<$end[$i] && $end1>$start[$i]) ||
            ($chr2==$chr[$i] && $start2<$end[$i] && $end2>$start[$i])) {
#	    print "MATCH!\n";
	    $output = ($i+1) . "\t$_";
	    print OUT $output;
	}
    }
}

rename $OUT_TEMP, $OUT;
close OUT;
close IN;
