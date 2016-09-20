use Symbol;

if ($#ARGV+1 != 1) {
    print "Usage: vcf2wig input.vcf\n";
    exit;
}

%sample_num = ();
%sample_out = ();
%sample_wchr = ();
%sample_wpos = ();

$ns = 0;

open(IN,"<$ARGV[0]");
while(<IN>) {
    next if m/^GENOME/;   # skip header lines

    @col = split(/\t/);
    # col[0] is "build"
    $chr = $col[1];
    $chr="23" if $chr eq "X";
    $chr="24" if $chr eq "Y";

    $pos = $col[2];
    $sample = $col[3];
    $cov = 0;
    for ($i=4;$i<=$#col;$i++) { $cov=1 if $col[$i]>0; }

    $num = $sample_num{$sample};
    if ($num==null) {     # new sample
	$num = ++$ns;
	$sample_num{$sample} = $num;
	$fh = gensym;
	open($fh,">$sample.wig");
	$sample_out{$sample} = $fh; 
	print $fh "track type=wiggle_0 name=Coverage viewLimits=0:1\n";
	$wchr = $chr;
	$wpos = $pos;
    }
    $out = $sample_out{$sample};
    $wchr = $sample_wchr{$sample};
    $wpos = $sample_wpos{$sample};

    if ($chr!=$wchr || $pos!=$wpos) {
	print $out "fixedStep chrom=chr$chr start=$pos step=1\n";
	$sample_wchr{$sample} = $chr;
    }

    print $out "$cov\n";
    $pos++;
    $sample_wpos{$sample} = $pos;
}

while( ($sample, $out) = each %sample_out) {
    close $out;
}

