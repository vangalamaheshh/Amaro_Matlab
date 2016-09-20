function [READN,READT,gcDistribN,gcDistribT]=normalize_read_weight_gc( READN, READT, referenceLane )

laneListN = unique(READN.lane);
laneListT = unique(READT.lane);

binGC = 0:1:100;
gcDistribN = zeros(length(binGC),length(laneListN));
gcDistribT = zeros(length(binGC),length(laneListT));

%---  Check that each read pair has G+C content  ---%
if ~isfield(READN,'gc')
    READN=calc_fragment_gc_sorted( READN );
end

if ~isfield(READT,'gc')
    READT=calc_fragment_gc_sorted( READT );
end


%---  Calculate G+C distribution in 1% windows  ---%
for i=1:length(laneListN)
    laneN = laneListN(i);
    idxLane = find(READN.lane==laneN);
    [gcDistribN(:,i),binGCN] = histc( READN.gc(idxLane), binGC );
    idxGCN(idxLane) = binGCN;
end
READN.idxGC=reshape(idxGCN,length(idxGCN),1);

for i=1:length(laneListT)
    laneT = laneListT(i);
    idxLane = find(READT.lane==laneT);
    [gcDistribT(:,i),binGCT] = histc( READT.gc(idxLane), binGC );
    idxGCT(idxLane) = binGCT;
end
READT.idxGC=reshape(idxGCT,length(idxGCT),1);


%---  Choose lane with median counts as the reference lane  ---%
if ~exist( 'referenceLane','var')
    [laneCounts,sortOrder] = sort( sum(gcDistribN) );
    referenceLane = sortOrder(floor(length(sortOrder)/2))
end


%---  Only correct G+C windows with at least 200 counts in normal lane  ---%
% GADDY suggests regularizing with pseudocount instead

minBinCount = 200;
[READN,gcWeightN] = calc_read_gc_weights( READN, gcDistribN, gcDistribN(:,referenceLane), minBinCount );

[READT,gcWeightT] = calc_read_gc_weights( READT, gcDistribT, gcDistribN(:,referenceLane), minBinCount );

