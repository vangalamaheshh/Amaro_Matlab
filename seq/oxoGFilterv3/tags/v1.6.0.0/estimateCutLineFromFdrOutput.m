function [x,y,minPoxCut, maxPoxNotCut] = estimateCutLineFromFdrOutput(fdrStruct, acs, PoxoG)
%
% x is a foxog value, y is the ac.
% TODO: Finish documentation.
%

pox = fdrStruct.pox(1:length(fdrStruct.cut));

[minPoxCut, minPoxCut_i] = min(pox(fdrStruct.cut));
[maxPoxNotCut, maxPoxNotCut_i] = max(pox(~fdrStruct.cut));

if isempty(maxPoxNotCut)
    maxPoxNotCut = 0;
end

if isempty(minPoxCut)
    minPoxCut = 1;
end

x = zeros(length(acs),1);
y = acs;
for i = 1:length(acs)
    ac = acs(i);
    
    acVals = [0:ac];
    tmpPox = binocdf(acVals, ac, PoxoG);
    acVal_i = find((tmpPox < minPoxCut) & (tmpPox > maxPoxNotCut),1,'last');
    if isempty(acVal_i)
       hi = find(tmpPox > maxPoxNotCut, 1);
       lo = find(tmpPox < minPoxCut, 1, 'last');
       if isempty(hi)
           hiVal = 0;
       else
           hiVal = acVals(hi);
       end
       if isempty(lo)
           loVal = 0;
       else
           loVal = acVals(lo);
       end
       x(i) = ( hiVal + loVal)/2;
    else
       x(i) = acVals(acVal_i);
    end
end
