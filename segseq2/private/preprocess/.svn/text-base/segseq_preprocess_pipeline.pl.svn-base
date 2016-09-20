#!/usr/local/bin/perl -w
#
#  FILE: segseq_preprocess_pipeline.pl

#  KNOWN ISSUES:
#  Path to Matlab code in splitFlowcells()

use strict;
use POSIX;

# Convert to getopt for optional parameters, e.g. qualCutoff
if ( @ARGV < 3 )
{
    die "Usage: perl segseq_preprocess_pipeline.pl <Normal read dir> <Tumor read dir> <Sample name> <Qual cutoff> <Exe dir> <Output dir>\n";
}


MAIN:
{
    my $normalDir = $ARGV[0] . "/";
    my $tumorDir = $ARGV[1] . "/";
    my $sampleName = $ARGV[2];
    my $qualCutoff = $ARGV[3];
    my $exeDir = $ARGV[4];
    my $baseDir = $ARGV[5];

    $baseDir = "/xchip/tcga_scratch2/ng/" . $sampleName unless ( defined $baseDir );
    my $outDir = $baseDir . "/dat/";
    my $matDir = $baseDir . "/matfiles/";

    $exeDir = "~/CancerGenomeAnalysis/trunk/matlab/segseq2/private/preprocess/" unless ( defined $exeDir );
    $qualCutoff = 20 unless ( defined $qualCutoff );
    
    # Enable Matlab compile libraries
    print STDERR "Making directory $outDir\n";
    `mkdir $baseDir`;
    `mkdir $outDir`;
    `mkdir $outDir/scripts`;
    `mkdir $baseDir/matfiles`;

    #---  1) CONVERT text files to sorted Matlab structures  ---#
    opendir( DATDIR, $outDir ) or die "ERROR: Can't open $outDir\n";
    my @datfiles = grep { $_ =~ /mat$/ } readdir ( DATDIR );
    
    if ( scalar @datfiles < 1 )
    {
	print STDERR "START conversion to Matlab structures\n";
	makeLSFscripts( $normalDir, $outDir, $exeDir, $qualCutoff, "normal" );
	makeLSFscripts( $tumorDir, $outDir, $exeDir, $qualCutoff, "tumor" );

	my $bjob = `bjobs`;
	my @jobList = split( "\n", $bjob );
	my $time = 0;
	while( scalar @jobList > 1 )
	{
	    print STDERR $time, " min\t", scalar( @jobList ), " LSF jobs remaining\n";
	    sleep( 60 );
	    $time += 1;
	    
	    $bjob = `bjobs`;
	    @jobList = split( "\n", $bjob );
	}
	print STDERR "FINISHED conversion to Matlab structures\n";
    }

    #---  2) Collect normal and tumor reads  ---#
    # FORK the Matlab process to compile reads
    my $matfileAll = $matDir ."/" . $sampleName . "_tumor_normal_aligned_paired_reads_qual" . $qualCutoff . ".mat";
 
#    if ( -x $matfileAll )
    {
        if ( ! defined( my $childID = fork() ) )
        {
   	    die "ERROR: Can't fork a sub-process in collectReads()\n";
        }
        elsif( $childID == 0 )
        {
	    print STDERR "Collecting normal reads\n";
   	    collectReads( $sampleName, $matDir, $outDir, $qualCutoff, "normal" );
	    print STDERR "Collecting tumor reads\n";
	    collectReads( $sampleName, $matDir, $outDir, $qualCutoff, "tumor" );

	    mergeTumorNormalReads( $matDir, $sampleName, $qualCutoff );

        }
        else
        {
	    waitpid( $childID, 0 );
	    print STDERR "Completed read collection\n";
        }
    }
    #---  3) Correct G+C  ---#
#    if (0 )
#    {
#     correctGC( $normalDir, $matDir, $sampleName, $qualCutoff );
#    }

    #---  4) Split flowcells  ---#
   splitFlowcells( $normalDir, $tumorDir, $matDir, $sampleName, $qualCutoff );

}

sub collectReads
{
    my ( $sampleName,
	 $matDir, 
	 $outDir,
	 $qualCutoff,
	 $type ) = @_;
    
    my $varName = "READT";
    $varName = "READN" if ( $type eq "normal" );

    my $mFile = ">" . $matDir . "segseq_preprocess_1_collect_reads_" . $type . ".m";
    open( OUT, $mFile ) or die "ERROR: Can't write to $mFile\n";

    print OUT "sampleName = \'", $sampleName, "\';\n";
    print OUT "type = \'", $type, "\';\n";
    print OUT "matfile = [ \'$matDir/\' sampleName \'_\' type \'_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n\n";
    print OUT "chrList = 1:23;\n\n";

    print OUT $varName . ".chr = [];\n";
    print OUT $varName . ".pos1 = [];\n";
    print OUT $varName . ".pos2 = [];\n";
    print OUT $varName . ".lane = [];\n\n";

    print OUT "for c=1:length(chrList)\n";
    print OUT "\tchr=chrList(c);\n";
    print OUT "\tfilename = [ \'", $outDir, "/" , $type, "_chr\' num2str(chr) '.txt.mat' ]\n";
    print OUT "\tR = load( filename );\n";

    print OUT "\t", $varName . ".chr = [ " . $varName . ".chr; R.P.chr ];\n";
    print OUT "\t", $varName . ".pos1 = [ " . $varName . ".pos1; R.P.pos1 ];\n";
    print OUT "\t", $varName . ".pos2 = [ " . $varName . ".pos2; R.P.pos2 ];\n";
    print OUT "\t", $varName . ".lane = [ " . $varName . ".lane; R.P.lane ];\n";

    print OUT "end\n\n";
    print OUT "save(matfile, \'", $varName, "\',\'-v7.3\');\n";    

    close( OUT );
    $mFile = substr( $mFile, 1 );

    `matlab -nodisplay < $mFile &`;

}


sub correctGC
{
    my ( $normalDir, 
	 $matDir,
	 $sampleName,
	 $qualCutoff ) = @_;


    my $mFile = ">" . $matDir . "/segseq_preprocess_3A_normalize_reads.m";
    open( OUT, $mFile ) or die "ERROR: Can't write to $mFile\n";

    print OUT "addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq/\n";
    print OUT "sampleName = \'", $sampleName, "\';\n";
    print OUT "matfileTN = [ \'$matDir/\' sampleName \'_tumor_normal_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n\n";
    print OUT "load(matfileTN);\n";

    print OUT "laneListN = unique(READN.lane);\n";
    print OUT "[numReads,binReads]=histc(READN.lane, laneListN);\n\n";
    print OUT "profile on\n";
    print OUT "[READN,READT,gcDistribN,gcDistribT] = normalize_read_weight_gc( READN, READT );\n\n";
    print OUT "profile off\n";

    print OUT "save(matfileTN, \'READN\',\'READT\',\'-v7.3\');\n\n";
    print OUT "profile report\n";

    close( OUT );
    $mFile = substr( $mFile, 1 );

    `matlab -nodesktop < $mFile &`;

}


sub makeLSFscripts
{
    my ( $dir,
	 $outDir,
	 $exeDir,
	 $qualCutoff,
	 $type ) = @_;
    opendir( DIR, $dir ) or die "ERROR: Can't open $dir\n";

    my $exe = $exeDir . "/convert_aligned_pairs_to_mat";


    my @files = grep { $_ =~ /chr/ } readdir( DIR );
    for my $f ( @files )
    {
	my $fullName = $dir . $f;
	my $matName = $outDir . $type . "_" . $f . ".mat";
	my $shFile = ">" . $outDir . "scripts/" . $type . "_" . $f . ".sh";
	open( OUT, $shFile ) or die "ERROR: Can't write $shFile\n";

	my $fileinfo = `ls -l $fullName`;
	my @filedetails = split( '\s', $fileinfo );

	my $memlimit = floor( $filedetails[4] / 1e9 );
	$memlimit = $memlimit * 4000;

	print OUT "LD_LIBRARY_PATH=/broad/tools/apps/matlab75/sys/os/glnxa64:/broad/tools/apps/matlab75/bin/glnxa64:/util/gcc-4.1.1/lib:\$LD_LIBRARY_PATH\n";
	print OUT "export LD_LIBRARY_PATH\n\n";

        print OUT "bsub -q priority -P tcga -R \"rusage[mem=$memlimit]\" -J $f -e err.txt -o log.txt $exe $fullName $matName $qualCutoff\n";


	close( OUT );
	$shFile = substr( $shFile, 1 );
	`chmod +x $shFile`;
	`$shFile`;
    }
}

sub mergeTumorNormalReads
{
    my ( $matDir,
	 $sampleName,
	 $qualCutoff ) = @_;

    my $matfileN = $matDir ."/" . $sampleName . "_normal_aligned_paired_reads_qual" . $qualCutoff . ".mat";
    my $matfileT = $matDir ."/" . $sampleName . "_tumor_aligned_paired_reads_qual" . $qualCutoff . ".mat";
    my $matfileAll = $matDir ."/" . $sampleName . "_tumor_normal_aligned_paired_reads_qual" . $qualCutoff . ".mat";

    
    my $mFile = ">" . $matDir . "segseq_preprocess_2_merge_tumor_normal" . ".m";
    open( OUT, $mFile ) or die "ERROR: Can't write to $mFile\n";

    print OUT "sampleName = \'$sampleName\'\n";
    print OUT "matfileN = [ \'$matDir/\' sampleName \'_normal_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n";
    print OUT "matfileT = [ \'$matDir/\' sampleName \'_tumor_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n\n";
    print OUT "matfileTN = [ \'$matDir/\' sampleName \'_tumor_normal_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n\n";

    print OUT "load(matfileN);\n";
    print OUT "load(matfileT);\n";
    print OUT "save(matfileTN,\'READN\',\'READT\',\'-v7.3\');\n";
    close( OUT );

    $mFile = substr( $mFile, 1 );
    `matlab -nodisplay < $mFile`;

}


sub splitFlowcells
{
    my ( $normalDir, 
	 $tumorDir, 
	 $matDir,
	 $sampleName,
	 $qualCutoff ) = @_;


    my $laneInfoN = $normalDir . "/lanelist.txt";
    my $laneInfoT = $tumorDir . "/lanelist.txt";
    my $laneListN = ">" . $matDir . "/laneinfoN.txt";
    my $laneListT = ">" . $matDir . "/laneinfoT.txt";
    my $laneList = ">" . $matDir . "/laneinfo.txt";

    my %flowcellCountN;
    my %flowcellIndexN;

    my %flowcellCountT;
    my %flowcellIndexT;

    my %flowcellInfo;

    print STDERR "WARNING: Assumes 0-indexing from MakeA.java but 1-indexing from MakeLaneList.java\n";
    print STDERR "Change offset to 0 if this discrepancy is changed\n";
    my $offset = 1;
    # CHANGE OFFSET TO ZERO 
    #my $offset = 0;


    #---  NORMAL LANES  ---#
    my $fcNum = 0;
    open( INFO, $laneInfoN ) or die "ERROR: Can't open file $laneInfoN\n";
    while( my $line = <INFO> )
    {
	chomp $line;
	my @row = split( "\t", $line );
	my $currFC = substr( $row[1], 0, 5 );
	my ( $currLane ) = $row[1] =~ /\S+\.([1-8])/;

	$flowcellCountN{$currFC}++;
	if ( $flowcellCountN{$currFC} == 1 )
	{
	    $fcNum++;
	    $flowcellIndexN{$currFC} = $fcNum;
	}
	my $infoString = join( "\t", $row[0] + $offset, $currFC, $fcNum, $currLane, 0 );
	$flowcellInfo{$currFC}{$currLane} = $infoString;
    }
    close( INFO );

    #---  TUMOR LANES  ---#
    open( INFO, $laneInfoT ) or die "ERROR: Can't open file $laneInfoT\n";
    while( my $line = <INFO> )
    {
	chomp $line;
	my @row = split( "\t", $line );
	my $currFC = substr( $row[1], 0, 5 );
	my ( $currLane ) = $row[1] =~ /\S+\.([1-8])/;

	if ( defined $flowcellIndexN{$currFC} )
	{
	    $fcNum = $flowcellIndexN{$currFC};
	}
	else
	{
	    $fcNum++;
	}
	my $infoString = join( "\t", $row[0] + $offset, $currFC, $fcNum, $currLane, 1 );
	$flowcellInfo{$currFC}{$currLane} = $infoString;

    }
    close( INFO );

    open( LANES, $laneList ) or die "ERROR: Can't open $laneList\n";
    print LANES "Lane number\tFlowcell name\tFlowcell number\tLane in flowcell\tNormal(0) | Tumor(1)\n";
    for my $fc ( sort keys %flowcellInfo )
    {
	for my $lane ( sort { $a <=> $b } keys %{$flowcellInfo{$fc}} )
	{
	    print LANES $flowcellInfo{$fc}{$lane}, "\n";
	}
    }
    close( LANES );

    
    my $mFile = ">" . $matDir . "segseq_preprocess_3_split_flowcells.m";
    open( OUT, $mFile ) or die "ERROR: Can't write to $mFile\n";

    print OUT "addpath ~/CancerGenomeAnalysis/trunk/matlab/segseq/\n";
    print OUT "sampleName = \'$sampleName\'\n";
    print OUT "matfile = [ \'$matDir/\' sampleName \'_tumor_normal_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n\n";
    print OUT "R = load(matfile);\n";

    print OUT "laneinfo = \'$matDir/laneinfo.txt\';\n";
    print OUT "fidN = fopen( laneinfo );\n";
    print OUT "INFO = textscan( fidN, \'%u8%s%u8%u8%u8\', \'headerLines\', 1 );\n";
    print OUT "FC = [ INFO{1} INFO{3} INFO{4} INFO{5} ];\n";
    print OUT "fcString = INFO{2};\n";
    print OUT "fcList = unique( FC(:,2) );\n\n";

    print OUT "fcReads = zeros( length(fcList), 1 );\n";
    print OUT "for f=1:length(fcList)\n";
    print OUT "    idxFC = find( FC(:,2) == fcList(f) );\n";
    print OUT "    currFC = FC( idxFC, : );\n";
    print OUT "    laneListN = currFC( find( currFC(:,4) == 0 ), 1 );\n";
    print OUT "    laneListT = currFC( find( currFC(:,4) == 1 ), 1 );\n";
    print OUT "    for n = 1:length(laneListN)\n";
    print OUT "        fcReads(fcList(f)) = fcReads(fcList(f)) + length( find(R.READN.lane == laneListN(n) ) );\n";
    print OUT "    end\n";
    print OUT "end\n\n";

    # Find flowcell with the maximum number of reads
    print OUT "totalReads=fcReads\n";
    print OUT "[sortedReads,idxSorted] = sort(fcReads)\n";
    print OUT "idxMedian = idxSorted(ceil(length(idxSorted)/2));\n";
    print OUT "    idxFC = find( FC(:,2) == idxMedian );\n";
    print OUT "    currFC = FC( idxFC, : );\n";
    print OUT "    laneListN = currFC( find( currFC(:,4) == 0 ), 1 );\n";
    print OUT "    laneListT = currFC( find( currFC(:,4) == 1 ), 1 );\n";
    print OUT "    READN = filter_reads_by_lane_num( R.READN, laneListN );\n";
    print OUT "    READT = filter_reads_by_lane_num( R.READT, laneListT );\n\n";
    print OUT "load HG18_N36_D2_WINDOWS_100K;\n";
    print OUT "RATIOS=calc_ratios_from_reads(READN,READT,WINDOWS);\n";

    print OUT "    currFCstring = fcString{ idxFC(1) }\n";
    print OUT "    outmat = [ \'$matDir/\' sampleName \'_\' currFCstring \'_tumor_normal_aligned_paired_reads_qual\' num2str(" . $qualCutoff . ") \'.mat\' ]\n\n";

    print OUT " save( outmat, \'READN\', \'READT\', \'RATIOS\', \'-v7.3\');\n";
    close( OUT );

    $mFile = substr( $mFile, 1 );
    `matlab -nodisplay < $mFile &`;

}

