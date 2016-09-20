function RATIOS=calc_ratios_from_reads(READN,READT,WINDOWS);

% Count number of reads in genome windows
currChrStartN = 1;
currChrStartT = 1;

numChr = max(READN.chr);
aN = length(READN.chr);
aT = length(READT.chr);

chr = [];
windows = [];
ratios = [];

for c=1:numChr
    currWindows = WINDOWS.breaks(find(WINDOWS.chr==c));

    currChrPosN = sort( READN.pos( find(READN.chr==c) ) );
    currChrPosT = sort( READT.pos( find(READT.chr==c) ) );

    currCountN = histc(currChrPosN, currWindows );
    currCountN = currCountN(1:(end-1));    % Last bin are exact matches

    currCountT = histc(currChrPosT, currWindows );
    currCountT = currCountT(1:(end-1));    % Last bin are exact matches

    currChr = repmat(c,length(currCountN),1);
    currPos = currWindows(1:length(currCountN));

    idxFilter = find(currCountN>0);
    currRatios = ( currCountT(idxFilter) / aT ) ./ ( currCountN(idxFilter) / aN );
    currChr = currChr(idxFilter);
    currPos = currPos(idxFilter);

    
    if c == 1
        chr = currChr;
        windows = currPos;
        ratios = currRatios;
    else
        chr = [ chr; currChr ];
        windows = [ windows; currPos ];
        ratios = [ ratios; currRatios ];
    end
end

RATIOS.chr = chr;
RATIOS.windows = windows;
RATIOS.ratios = ratios;
