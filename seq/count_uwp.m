function nuwp = count_uwp(R,idx,chr,pos,params)
% count number of unique weird-pair destinations
% subfunction of somcall() and germcall()
% Mike Lawrence 2009-08

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'weirdpair_threshold',10000);
params = impose_default_value(params,'weirdpair_distinctness_threshold',2000);

wp = R(idx,10:12);                                                     % chr start strand
wp = wp(wp(:,1)>0,:);                                                  % remove unpaired and unmapped-pair
nw = find(wp(:,1)==chr & abs(wp(:,2)-pos)<params.weirdpair_threshold); % remove non-weird
wp = wp(setdiff(1:size(wp,1),nw),:);
wp(:,2) = round(wp(:,2) / params.weirdpair_distinctness_threshold);    % collapse to same destination
nuwp = size(unique(wp,'rows'),1);
