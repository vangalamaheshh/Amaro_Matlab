$usage = 'bsub "cat results_sorted.txt | perl joindump.pl [<num_cols_to_match,default=2>] > results_joined.txt"';

if ($#ARGV+1 == 0) {
    $NTM = 2;
} elsif ($#ARGV+1 == 1) {
    $NTM = $ARGV[0];
} else {
    die($usage);
}

@fields1 = (); @fields2 = ();
while (<STDIN>) {
  chomp;
  @fields2 = split(/\t/);
  $n = $#fields1;
  if ($n>=$NTM && $#fields2==$n) {
    $is_match = 1;
    for ($i=0;$i<$NTM;$i++) {
	if ($fields1[$i] ne $fields2[$i]) {
	    $is_match = 0;
	    last;
	}
    }
    if ($is_match) {
	for(;$i<=$n;$i++) {
	    if ($fields1[$i]==-1) { $fields1[$i]=$fields2[$i]; }
	}
	@fields2 = ();
    }
  }
  for ($i=0;$i<=$n;$i++) {
     print $fields1[$i];
     if ($i<$n) { print "\t" } else { print "\n" }
  }
  @fields1 = @fields2;
}
