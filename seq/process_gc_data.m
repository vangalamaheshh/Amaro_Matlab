function T = process_gc_data(iszdir)
% process_gc_data(iszdir)
%
% Process GCContentByLane data
%
% Mike Lawrence 2009-07-30

T = [];
T.tally_desc = '1=read only(actual) 2=whole insert(ref)';

fprintf('Chromosomes: ');for c=1:24, fprintf('%d/24 ',c);
  tmp = load_matrix([iszdir '/chr' num2str(c) '.gc']);
  tmp = cat(3,tmp(1:101,3:end),tmp(102:202,3:end));
  if c==1, T.raw = tmp; else T.raw(:,:,:,c) = tmp; end
end, fprintf('\n');
T.raw_desc = 'gc% x lane x tally x chr';

T.dat = sum(T.raw,4);
T.dat_desc = 'gc% x lane x tally';

T.nlanes = size(T.dat,2);

T.cov = squeeze(sum(T.dat,1));
T.good = find(T.cov(:,1)>=0.1*median(T.cov(:,1)) ...
        & T.cov(:,2)>=0.1*median(T.cov(:,2)));

% normalize each lane
T.datnorm = T.dat ./ cat(3,repmat(T.cov(:,1)',101,1),repmat(T.cov(:,2)',101,1));
T.datnorm(:,setdiff(1:T.nlanes,T.good),:) = nan;
