#!/usr/bin/perl
# create session XML and zipped seg files

use strict;
use warnings;

my $usage = "Usage: $0 seg_in_dir seg_out_dir [session_out_dir seg_out_URL]";

my $nargs = @ARGV;

if ($nargs < 2) {
    die("$usage");
}


## mandatory arguments

# IGV segment file (input) directory
my $segdir = shift;  #!"/xchip/gistic/tcgascape/tcgascape_120416/igvfiles";

# zipped segment file output directory (currently served as an URL)
my $datapath = shift; #! "/xchip/gistic/tcgascape/igvfiles";


## optional arguments

my $session_outdir = "";
my $servurl = "";
if ($nargs > 3) {
    # session xml output directory
    $session_outdir = shift; #!"/xchip/gistic/tcgascape/tcgascape_120416/session_xml";

    # URL for accessing $datapath via http
    $servurl = shift; #!"http://www.broadinstitute.org/igvdata/tcga/tcgascape";
}
# iterate across all files, generate xml and zip for *.seg
opendir(DIR, $segdir) or die $!;
while ($_ = readdir(DIR)) {
    if( m/^(.+)\.seg$/ ) {
	my $ctype = lc($1);
	$ctype =~ tr/ /_/;
	if( $1 ne $ctype ) {
	    system("mv '$1.seg' $ctype.seg");
	}
	system("gzip -c $segdir/$ctype.seg > $datapath/$ctype.seg.gz");
	# write out session files, if directory was supplied
	if ($session_outdir) {
	    write_xml($ctype);
	}
    }
}
closedir(DIR);

# subroutine to write session file for a given cancer type
sub write_xml {
    my $ctype = shift;
    open(XFILE,">$session_outdir/${ctype}_session.xml") or die;
    print XFILE <<END;
<?xml version="1.0" encoding="UTF-8"?>
<Global genome="hg19" version="1">
    <Resources>
        <Resource path="$servurl/$ctype.seg.gz"/>
        <Resource path="$servurl/sample_info.txt"/>
    </Resources>
</Global>
END
    close XFILE;
}
