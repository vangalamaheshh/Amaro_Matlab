function T = dRanger_examine_twinned_breakpoints(E,J)
% T = dRanger_examine_twinned_breakpoints(E,J)
% E and J come from dRanger_find_chains(x,P)

nxi = size(J,1)/2;
if nxi~=round(nxi), error('J should have even number of rows, same number of columns'); end

% look at all inter-pair adjancencies, see whether there is overlap
J2 = J; for j=1:nxi, J2(j,j+nxi)=0; J2(j+nxi,j)=0; end   % remove within-rearrangement adjacencies
for j=1:nxi*2, J2(j,1:j) = 0; end  % remove diagonal and below
[a b] = find(J2);
T = sortrows([E(a,4:6) E(b,5:6)]);
T(:,end+1) = T(:,4)-T(:,2);
T(:,end+1) = 0;
T(T(:,3)==0 & T(:,5)==1 & T(:,4)>T(:,2),end) = 1;
T(T(:,3)==1 & T(:,5)==0 & T(:,2)>T(:,4),end) = 1;
T(T(:,3)==T(:,5),end) = nan;
