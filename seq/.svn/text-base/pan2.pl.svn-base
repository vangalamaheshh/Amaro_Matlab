#!/usr/bin/env perl

die "usage:   pan2 outstem\n" if ($#ARGV+1 != 1);
$outstem = $ARGV[0];

# Firehose API Guide
# http://iwww.broadinstitute.org/cancer/cga/wiki/index.php/Firehose_API_Guide

$loginfile = "/home/unix/lawrence/script/fhlogin";
$curl = "curl -K $loginfile -sL";

# columns to scrape for each individual
@icolnames = ('individual_id','bsp_participant_id',
              'maf_file_capture','indel_maf_file_capture','somatic_mutation_coverage_capture',
  	      'maf_file_wgs','indel_maf_file_wgs','somatic_mutation_coverage_wgs',
 	      'maf_file_oxoG_capture','maf_file_capture2','somatic_mutation_coverage_capture2');
$numicols = $#icolnames+1;

# columns to scrape for each sample
@scolnames = ('sample_id','bam_file_capture','bam_file_wgs','clean_bam_file_capture',
              'clean_bam_file_wgs');
$numscols = $#scolnames+1;

# open output files and write header lines

$ioutfile = "$outstem.individuals.txt";
$soutfile = "$outstem.samples.txt";

print "Output files\n\tIndividuals: $ioutfile\n\tSamples: $soutfile\n";
open(iOUT, ">$ioutfile");
print iOUT "workspace\tbuild\tindividual_set";
for($ci=0;$ci<$numicols;$ci++) {
    print iOUT "\t$icolnames[$ci]";
}
print iOUT "\n";

print sOUT "\n";
open(sOUT, ">$soutfile");
print sOUT "workspace\tbuild\tindividual_set";
for($ci=0;$ci<$numscols;$ci++) {
    print sOUT "\t$scolnames[$ci]";
}
print sOUT "\n";

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
    next if ($wno<413);
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

    # iterate through individual sets

    $ict = $#isets+1;
    $ino = 0;
    foreach $iset (@isets) {
	$ino++;
	chomp $iset;

	printf "(workspace %7s)  %-47s (iset %5s)  %-47s",
	     sprintf("%d/%d",$wno,$wct),$workspace,sprintf("%d/%d",$ino,$ict),$iset;

	# get list of individuals in this individual set

	$cmd = "$curl 'http://firehose:8080/cga/ws/entity/getAnnotations/Individual?".
	    "entityNames=$iset&filterSetType=Individual_Set&workspaceName=$workspace'";
	$result = `$cmd`;
	@lines = split('^',$result);
	if ($lines[0] =~ /<html>/) {
	    print "Error getting individual list\n";
	    next;
	}

	# find relevant columns
	@columns = split('\t',$lines[0]);
	shift @lines;  # remove header line
	my @colidx;
	for($ci=0;$ci<$numicols;$ci++) {
	    $colidx[$ci] = find_string($icolnames[$ci],\@columns);
	}
    
	# output <workspace_id> <build> <individual_set_id> <individual_id> (mafs... wigs...)
	$iict = 0;
	foreach $line (@lines) {
	    @fields = split('\t',$line);
	    print iOUT "$workspace\t$build\t$iset";
	    for($ci=0;$ci<$numicols;$ci++) {
		$cell = '.';
		if ($colidx[$ci]!=-1) {
		    $cell = $fields[$colidx[$ci]];
		}
		chomp $cell;
		print iOUT "\t$cell";
	    }
	    print iOUT "\n";
	    $iict++;

	}  # next individual

	printf("[%5d]\n",$iict);   # display how many individuals

	# get list of samples in this individual set

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
	for($ci=0;$ci<$numscols;$ci++) {
	    $colidx[$ci] = find_string($scolnames[$ci],\@columns);
	}

        # output <workspace_id> <build> <individual_set_id> <sample_id> (bams...)
        $iict = 0;
        foreach $line (@lines) {
            @fields = split('\t',$line);
            print sOUT "$workspace\t$build\t$iset";
            for($ci=0;$ci<$numscols;$ci++) {
                $cell = '.';
                if ($colidx[$ci]!=-1) {
                    $cell = $fields[$colidx[$ci]];
                }
                chomp $cell;
                print sOUT "\t$cell";
            }
            print sOUT "\n";
	    
        }  # next sample

    }  # next individual_set
}  # next workspace

close iOUT;
close sOUT;




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



