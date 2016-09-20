# usage: bsub "cat results_sorted.txt | perl joinweird.pl > results_joined.txt"
$pre1 = 'xxx'; $post1 = 'xxx';
$n = 0; $n_pairs = 0; $n_unpaired = 0;
while (<>) {
  $n++;
  $result = (($pre2,$post2) = m/^(.*)\t([^\t]*\t[^\t]*)$/);
  if ($result==2) {
    if ($pre1 eq $pre2) {
      # join "post" components
      $result1 = (($ed11,$ed21) = ($post1 =~ m/^([-\d]*)\t([-\d]*)$/)); 
      $result2 = (($ed12,$ed22) = ($post2 =~ m/^([-\d]*)\t([-\d]*)$/));
      if ($result1==2 && $result2==2) {
        $ed1 = $ed11; if ($ed1==-1) { $ed1 = $ed12; }
        $ed2 = $ed21; if ($ed2==-1) { $ed2 = $ed22; }
        $post = "$ed1\t$ed2";
        print $pre1 . "\t" . $post . "\n";
      }
      else {
	  # print "ERROR\n";
      }
      $pre1 = 'xxx'; $post1 = 'xxx';
      $n_pairs++;
    } else {
      if ($pre1 ne "xxx") { $n_unpaired++; } #print "UNPAIRED\n" }
      $pre1 = $pre2; $post1 = $post2;
    }
  }
}
