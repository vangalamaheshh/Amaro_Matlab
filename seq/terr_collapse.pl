if ($#ARGV+1 != 1) {
    print "Usage: terr_collapse rulefile\n";
    exit;
}

$RULEFILE = $ARGV[0];
open (RULEFILE, "<$RULEFILE");
$max = 0;
while(<RULEFILE>) {
    chomp;
    $rule[++$i]=$_;
    $max = $_ if $_>$max;
}
close RULEFILE;

while(<STDIN>) {
    @col = split(/\t/);
    $#tot = -1;
    print "$col[0]\t$col[1]\t$col[2]\t$col[3]\t$col[4]";
    for ($i=5;$i<=$#col;$i++) {
	$old = $i-4;
	$new = $rule[$old];
	$tot[$new] += $col[$i];
    }
    for ($i=1;$i<=$max;$i++) {
	print "\t$tot[$i]";
    }
    print "\n";
}




