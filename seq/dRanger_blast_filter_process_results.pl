$INTABLE = $ARGV[0];
$HITSTEM = $ARGV[1];
$OUTTABLE = $ARGV[2];

$| = 1;  # piping-hot output

# load input table

if (!(-e $INTABLE)) { print "Input table $INTABLE not found!\n"; exit; }
open(IN,"<$INTABLE");
$j=0;
while (<IN>) {
    $result = (($jump[$j],$file[$j],$chr1[$j],$mn1[$j],$mx1[$j],
                                    $chr2[$j],$mn2[$j],$mx2[$j]) = 
	       m/^([-.\d]+)\t([-.\d]+)\t([-.\d]+)\t([-.\d]+)\t([-.\d]+)\t([-.\d]+)\t([-.\d]+)\t([-.\d]+)$/);
    if (!$result) { print "Invalid line in INTABLE:\n$_\n"; exit; }
    $j++;
}
close IN;
$nj = $j;

# check input table

if (!$nj) { print "Empty input table!\n"; exit; }
$oldf = 1;
for($j=0;$j<$nj;$j++) {
    if ($jump[$j] != $j+1) { print "jump $j should be number ".($j+1).", not $jump[$j]\n"; exit; }
    if ($file[$j] != $oldf) {
	if ($file[$j] == $oldf+1) { $oldf = $file[$j]; }
	else { print "files have improper order in input table\n"; exit; }
    }
}

# process files

@hits1 = ((0)x$nj);
@hits2 = ((0)x$nj);
@self1 = ((0)x$nj);
@self2 = ((0)x$nj);
@other1 = ((0)x$nj);
@other2 = ((0)x$nj);

for ($currfile=1;$currfile<=$file[$nj-1];$currfile++) {
    $currfilename = $HITSTEM . $currfile;
    print "Analyzing file $currfile of $file[$nj-1]\n";
    if (!(-e $currfilename)) { print "Hit file $currfilename not found!\n"; exit; }
    open(IN,"<$currfilename");
    while (<IN>) {
	if (m/^#/) { next; }   # skip header lines
	chomp;
        $result = (($j,$e,$c,$st,$en) = 
		   m/^jump(\d+)\.end(\d)\tchr([\dXY]+).*\t(\d+)\t(\d+)\t[^\t]+\t[^\t]+$/);
        if (!$result) { print "REJECTED LINE: $_\n"; next; }
	if ($en<$st) { $tmp=$st; $st=$en; $en=$tmp; }  # make sure $st<$en
        if ($c eq "X") { $c = 23; } elsif ($c eq "Y") { $c = 24; }
	if ($j<1 || $j>$nj) { print "jump $j out of range!\n"; exit; }
	$j--; # switch to zero-based jump index
	if ($e==1) {
            $hits1[$j]++;
	    if ($c==$chr1[$j] && $st<=$mx1[$j] && $en>=$mn1[$j]) { $self1[$j]++; }
	    if ($c==$chr2[$j] && $st<=$mx2[$j] && $en>=$mn2[$j]) { $other1[$j]++; }
	} elsif ($e==2) {
	    $hits2[$j]++;
            if ($c==$chr2[$j] && $st<=$mx2[$j] && $en>=$mn2[$j]) { $self2[$j]++; }
            if ($c==$chr1[$j] && $st<=$mx1[$j] && $en>=$mn1[$j]) { $other2[$j]++; }
	} else {
	    print "invalid end $e\n"; exit;
	}
    }
    close(IN);
}

# write output table

print "Generating output table\n";
$OUTTEMP = $OUTTABLE . "_partial";
open(OUT,">$OUTTEMP");
for ($j=0;$j<$nj;$j++) {
  print OUT "$jump[$j]\t$hits1[$j]\t$self1[$j]\t$other1[$j]\t".
      "$hits2[$j]\t$self2[$j]\t$other2[$j]\n";
}
close(OUT);
rename $OUTTEMP, $OUTTABLE;
print "Done\n";


