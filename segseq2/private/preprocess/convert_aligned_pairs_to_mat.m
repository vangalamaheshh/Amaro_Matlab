function convert_aligned_pairs_to_mat(in_fname, out_fname, qualCutoff)

% Sort by lane ID, read ID, position
%<LaneID> <ReadID> <Chr> <Start> <End> <Strand> <MappingQuality> <WhichPairmate>

f=fopen(in_fname);

colLaneID = 1;
colReadID = 2;
colChr = 3;
colLeft = 4;
colRight = 5;
colQual = 7;

qualCutoff=str2num(qualCutoff)


C=dlmread(in_fname);

idxFilter = find(C(:,colQual)>= qualCutoff);
C=C(idxFilter,:);

C=sortrows(C,[colLaneID colReadID colChr colLeft]);

% MATCH READ PAIRS in sequential order
idxRead2 = find( diff(C(:,2))==0 ) + 1;

numReads = length(idxRead2);

% Reshape vector
P.chr = uint8(C(idxRead2, colChr));
P.pos1 = C(idxRead2-1, colLeft);
P.pos2 = C(idxRead2, colRight);
P.lane = single(C(idxRead2, colLaneID));

save(out_fname,'P','-v7.3');
