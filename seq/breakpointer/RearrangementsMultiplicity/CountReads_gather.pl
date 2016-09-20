#!/util/bin/perl -w
#
# Yotam Drier, yotamd@gmail.com
#
use strict;

die "Usage: CountReads_gather <sample.id> <libdir> <resultsdir> [<keeptempfiles>]\n" if ($#ARGV+1 < 3);
my $sample = $ARGV[0];
my $libdir = $ARGV[1];
my $resultsdir = $ARGV[2];
my $keeptempfiles = 0; 
if ($#ARGV+1 >= 4) {
    $keeptempfiles = ($ARGV[3] eq "1" || $ARGV[3] eq "true");
}

my ($start, $finish, $i, $dir, $line, $normsamp);
my @t;
my @dirs = `ls -dX1 $resultsdir/scatter.*/`;
my $results_manysamp="support.onenormal.txt";
my $results_onesamp=$sample.".support.txt";

if ($sample ne "NO_DATA") {
  chomp($dirs[0]);
  if (-e $dirs[0]."/$results_onesamp") {
    open(ALLBKPTS_FILE, "> $resultsdir/$results_onesamp");
    my $first = 1;
    foreach $dir (@dirs) {
      chomp($dir);	
      open(INFILE, "< $dir/$results_onesamp");
      if (!$first) {
	$line = <INFILE>;
      }
      $first = 0;
      foreach $line (<INFILE>) {		  
	print ALLBKPTS_FILE $line;
      }
      close(INFILE);
      if ($keeptempfiles==0) {
	system("rm -rf $dir/*/");
	system("rm $dir/$results_onesamp");
      }
    }
    close(ALLBKPTS_FILE);
  } else {
    foreach $dir (@dirs) {
      chomp($dir);	
      my @samples = `ls -X1 $dir`;
      foreach $normsamp (@samples) {
	chomp($normsamp);	
	if (-e "$dir/$normsamp/$results_manysamp") {
	  if (!-e "$resultsdir/$normsamp") {
	    mkdir "$resultsdir/$normsamp";
	    if ($keeptempfiles) {
	      system("cp $dir/$normsamp/$results_manysamp $resultsdir/$normsamp");
	    } else {
	      system("mv $dir/$normsamp/$results_manysamp $resultsdir/$normsamp");
	    }
	  } else {
	    open(ALLBKPTS_FILE, ">> $resultsdir/$normsamp/$results_manysamp");
	    open(INFILE, "< $dir/$normsamp/$results_manysamp");
	    $line = <INFILE>;
	    foreach $line (<INFILE>) {		  
	      print ALLBKPTS_FILE $line;
	    }
	    close(INFILE);
	    close(ALLBKPTS_FILE);
	  }
	} else {
	  die "Missing input file $dir/$normsamp/$results_manysamp";
	}
      }
      system("rm -rf $dir/*/") if $keeptempfiles==0;	  
    }
  }
}
