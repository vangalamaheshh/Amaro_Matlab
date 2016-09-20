# usage: bsub "cat chr1.txt | perl midchunks.pl 400 > chr1_value_at_center_of_400bp_chunks.txt"
$chunksize = $ARGV[0];
if ($chunksize<2 || $chunksize>100000000) {
    die("sumchunks parameters should be an integer from 2-100,000,000");
}

$midchunk = int($chunksize/2);

$numlines = 0;
while (<STDIN>) {
  chomp;
  $numlines++;
  if ($numlines == $midchunk) {
      print "$_\n";
  }
  if ($numlines == $chunksize) {
     $numlines = 0;
  }
}
