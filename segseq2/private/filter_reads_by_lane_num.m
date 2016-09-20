function READS=filter_reads_by_lane_num( IN_READS, laneList )

READS.chr=[];
READS.pos=[];
READS.lane=[];
idxKeep = [];

for i=1:length(laneList)
    laneNum = laneList(i);
    fprintf(1, [ 'Lane ' num2str(laneNum) '..' ] );
    idxLane = find(IN_READS.lane==laneNum);
    idxKeep = [ idxKeep; idxLane ];
end
fprintf(1,['\nTotal reads: ' num2str(length(idxKeep)) '\n' ]);
READS.chr = IN_READS.chr(idxKeep);
READS.pos = IN_READS.pos1(idxKeep);
READS.lane = IN_READS.lane(idxKeep);

compoundpos = READS.pos + 1e9*double(READS.chr);
[sortpos,sortOrder] = sort(compoundpos);
READS.chr=READS.chr(sortOrder);
READS.pos=READS.pos(sortOrder);
READS.lane=READS.lane(sortOrder);



