function R = load_aligned_single_reads( fname, qualCutoff )
%------------------------------------------------------------------------%
%  FILE: load_aligned_single_reads.m                                     %
%------------------------------------------------------------------------%
%                  SOFTWARE COPYRIGHT NOTICE AGREEMENT                   %
%      This software and its documentation are copyright (2008) by the   %
%  Broad Institute/Massachusetts Institute of Technology.  All rights    %
%  are reserved.  This software is supplied without any warranty or      %
%  guaranteed support whatsoever. Neither the Broad Institute nor MIT    %
%  can be responsible for its use, misuse, or functionality.             %
%------------------------------------------------------------------------%
%  INPUT:  Text file with coordinates of aligned reads                   %
%          Column 1 - Chr, Column 2 - Alignment start, Column 3 - Strand %
%------------------------------------------------------------------------%
%  PARAMETER LIST                                                        %
%    fname           Filename of text file with aligned read coordinates %
%    WINDOWS         Start coordinates for windows of equal size in the  %
%                      alignable portion of the reference genome         %
%      .chr          Chromosome                                          %
%      .breaks       Start coordinates of windows                        %
%------------------------------------------------------------------------%
%  OUTPUT VARIABLES                                                      %
%    COUNTS          Matlab structure with counts of aligned reads       %
%      .chr          Chromosome                                          %
%      .breaks       Start coordinates of windows                        %
%                      NOTE: length(COUNTS.breaks) = WINDOWS.breaks - 23 %
%                      because last bin of histc equals 0                %
%      .counts       Counts of aligned reads                             %
%                                                                        %
%    READS           Matlab structure with aligned read coordinates      %
%      .chr          Chromosome                                          %
%      .pos          Start positions of aligned reads                    %
%------------------------------------------------------------------------%

% Load coordinates of aligned position files with these column headings
%<LaneID> <ReadID> <Chr> <Start> <End> <Strand> <MappingQuality> <WhichPairmate>

fid = fopen(char(fname));
if fid < 0
    msg = [ 'ERROR: File ' char(fname) ' does not exist'];
    error('load_aligned_reads:fileNotFound',msg );
end

C = textscan(fid,'%u8%f64%u8%f64%f64%u8%u8%u8');
fclose(fid);

% Assign contents to variables 
laneID = C{1};
readID = C{2};
chrUnsorted = C{3};
posStart = C{4};
posEnd = C{5};
strand = C{6};
mapQuality = C{7};
whichPair = C{8};

idxFilter = find(whichPair==1 & mapQuality>= qualCutoff);

laneNum = laneID(idxFilter);
chrUnsorted = chrUnsorted(idxFilter);
posStart = posStart(idxFilter);
posEnd = posEnd(idxFilter);
strand = strand(idxFilter);
mapQuality = mapQuality(idxFilter);

chrList = unique(chrUnsorted);

sortpos = zeros(length(posStart),1);
sortchr = zeros(length(posStart),1,'single');
sortstr = zeros(length(posStart),1,'single');
sortlane = zeros(length(posStart),1,'single');
clear C;

% SORT ALIGNED START POSITIONS
chrBreaks = ones(length(chrList),1);
currChrStart = 1;

numReads = length(posStart);

for c=1:length(chrList)
    chr=chrList(c);
    idxChr = find(chrUnsorted==chr);

    [currChrPos,currOrder] = sort( posStart( idxChr ) );
    currChrEnd = currChrStart + length(currChrPos) - 1;
    sortpos(currChrStart:currChrEnd) = currChrPos;

    currChrStrand = strand( idxChr );
    currChrStrand = currChrStrand(currOrder);
    sortstr(currChrStart:currChrEnd) = currChrStrand;

    currLane = laneNum( idxChr );
    currLane = currLane(currOrder);
    sortlane(currChrStart:currChrEnd) = currLane;

    sortchr(currChrStart:currChrEnd) = chr;
    currChrStart = currChrStart + length(currChrPos);  % Update current start position
end

% CONVERT READ NAMES TO LANES
R.chr = sortchr;
R.pos1 = sortpos;
R.strand = sortstr;
R.lane = sortlane;

