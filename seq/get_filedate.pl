#!/usr/bin/env perl

die "usage:   check_exist filenamelist resultsfile\n" if ($#ARGV+1 != 2);
$infile = $ARGV[0];
$outfile = $ARGV[1];

open(IN, "<$infile");
open(OUT, ">$outfile");

$tot = 0; $y = 0; $n = 0;
while (<IN>) {
    chomp $_;
    if (-e $_) {
	open(HANDLE,$_);
	my $date = localtime( (stat HANDLE)[9] );
	close HANDLE;
	print OUT "$date\n";
	$y++;
    } else {
	print OUT "0\n";
	$n++;
    }
    $tot++;
    if ($tot % 100 == 0) {
	print "$y exist  +  $n nonexist  =  $tot total\n";
    }
}


close OUT;
