#!/usr/bin/env perl

die "usage:   panbam output_filename\n" if ($#ARGV+1 != 1);
$outfile = $ARGV[0];

# Firehose API Guide
# http://iwww.broadinstitute.org/cancer/cga/wiki/index.php/Firehose_API_Guide

$loginfile = "/home/unix/lawrence/script/fhlogin";
$curl = "curl -K $loginfile -sL";

# columns to scrape for each sample

@colnames = ('sample_id','bam_file_capture','bam_file_wgs','clean_bam_file_capture',
	     'clean_bam_file_wgs');
$numcols = $#colnames+1;

# open output file and write header line

print "Output file = $outfile\n";
open(OUT, ">$outfile");
print OUT "workspace\tbuild\tindividual_set";
for($ci=0;$ci<$numcols;$ci++) {
    print OUT "\t$colnames[$ci]";
}
print OUT "\n";

# get list of workspaces

print "Getting list of workspaces\n";
$cmd = "$curl 'http://firehose:8080/cga/ws/list/Workspace'";
$result = `$cmd`;
@workspaces = split('^',$result);
if ($workspaces[0] =~ /<html>/) {
    die "Error getting workspace list\n";
}

# survey each workspace

$wct = $#workspaces+1;
$wno = 0;
foreach $workspace (@workspaces) {
    $wno++;
    chomp $workspace;

    # find out what build this workspace is
    $cmd = "$curl 'http://firehose:8080/cga/ws/entity/getAnnotations/Workspace?entityNames=$workspace&workspaceName=$workspace'";
    $result = `$cmd`;
    if (result =~ /<html>/) {
        print "Error getting workspace annotations for workspace $workspace\n";
        next;
    }
    if ($result =~ /Homo_sapiens_assembly18/ && $result =~ /Homo_sapiens_assembly19/) {
	$build = 'hg18/hg19';
    } elsif ($result =~ /Homo_sapiens_assembly18/) {
	$build = 'hg18';
    } elsif ($result =~ /Homo_sapiens_assembly19/) {
        $build = 'hg19';
    } else {
        $build = 'other';
    }

    # get list of individual_sets

    $cmd = "$curl 'http://firehose:8080/cga/ws/list/Individual_Set?workspaceName=$workspace'";
    $result = `$cmd`;
    @isets = split('^',$result);
    if ($isets[0] =~ /<html>/) {
	print "Error getting individual_set list for workspace $workspace\n";
	next;
    }

    # get list of samples in each individual set

    $ict = $#isets+1;
    $ino = 0;
    foreach $iset (@isets) {
	$ino++;
	chomp $iset;

	printf "(workspace %7s)  %-47s (iset %5s)  %-47s",
	     sprintf("%d/%d",$wno,$wct),$workspace,sprintf("%d/%d",$ino,$ict),$iset;


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
	my @colidx;
	for($ci=0;$ci<$numcols;$ci++) {
	    $colidx[$ci] = find_string($colnames[$ci],\@columns);
	}
    
	# output <workspace_id> <build> <individual_set_id> <sample_id> (bams...)
	$iict = 0;
	foreach $line (@lines) {
	    @fields = split('\t',$line);
	    print OUT "$workspace\t$build\t$iset";
	    for($ci=0;$ci<$numcols;$ci++) {
		$cell = '.';
		if ($colidx[$ci]!=-1) {
		    $cell = $fields[$colidx[$ci]];
		}
		chomp $cell;
		print OUT "\t$cell";
	    }
	    print OUT "\n";
	    $iict++;

	}  # next sample

	printf("[%5d]\n",$iict);   # how many samples

    }  # next individual_set
}  # next workspace

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



