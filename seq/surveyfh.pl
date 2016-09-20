#!/usr/bin/env perl

die "usage:   surveyfh workspace_id output_filename\n" if ($#ARGV+1 != 2);

$workspace = $ARGV[0];
print "Surveying workspace $workspace\n";

$outfile = $ARGV[1];
print "Output file = $outfile\n";
open(OUT, ">$outfile");

$loginfile = "/home/unix/lawrence/script/fhlogin";
$curl = "curl -K $loginfile -sL";

# get list of individual_sets

$cmd = "$curl 'http://firehose:8080/cga/ws/list/Individual_Set?workspaceName=$workspace'";
$result = `$cmd`;
@isets = split('^',$result);

if ($isets[0] =~ /<html>/) {
    print "Error getting individual_set list\n";
    exit(1);
}

# get list of samples in each 

foreach $iset (@isets) {
    chomp $iset;
    next if ($iset !~ /^PR/);     # only consider production sets

    # get sample list
    print "Surveying individual_set $iset\n";
    $cmd = "$curl 'http://firehose:8080/cga/ws/entity/getAnnotations/Sample?".
	"entityNames=$iset&filterSetType=Individual_Set&workspaceName=$workspace'";
    $result = `$cmd`;
    @lines = split('^',$result);
    if ($lines[0] =~ /<html>/) {
	print "Error getting sample list\n";
	next;
    }

    # find relevant columns
    @columns = split('\t',$lines[0]);
    shift @lines;  # remove header line
    $id_col = find_string('sample_id',\@columns);
    $cbam_col = find_string('bam_file_capture',\@columns);
    $wbam_col = find_string('bam_file_wgs',\@columns);
    if ($id_col==-1) {
	print "No sample_id column\n";
	next;
    }
    if ($cbam_col==-1 && $wbam_col==-1) {
        print "No bam_file columns\n";
	next;
    }
    
    # output <individual_set_id>  <sample_id>  <wgs/capture>  <bamfile>
    foreach $line (@lines) {
	@fields = split('\t',$line);
	chomp($id = $fields[$id_col]);
	chomp($cbam = $fields[$cbam_col]);
	chomp($wbam = $fields[$wbam_col]);
	next if ($id eq "");
	print OUT "$iset\t$id\tcapture\t$cbam\n" if ($cbam ne "");
        print OUT "$iset\t$id\twgs\t$wbam\n" if ($wbam ne "");
    }
}

close OUT;





sub find_string {
    my ($target,$ref_strings) = @_;
    @strings = @{$ref_strings};

    $target = lc($_[0]);
    $i=0;
    foreach $string (@strings) {
	return $i if (lc($string) eq $target);
	$i++;
    }
    return -1;
}



