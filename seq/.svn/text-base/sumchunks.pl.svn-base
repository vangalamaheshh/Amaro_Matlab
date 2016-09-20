# usage: bsub "cat chr1.cbb | perl sumchunks.pl 400 > chr1_coverage_by_chunks_of_400.txt"
$chunksize = $ARGV[0];
if ($chunksize<2 || $chunksize>100000000) {
    die("sumchunks parameters should be an integer from 2-100,000,000");
}

$tot = 0;
$numlines = 0;
while (<STDIN>) {
  chomp;
  $tot += $_;
  $numlines++;
  if ($numlines == $chunksize) {
      print "$tot\n";
      $tot = 0;
      $numlines = 0;
  }
}
