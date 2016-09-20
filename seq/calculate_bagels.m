function M = calculate_bagels(M,P)

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'bagel_size', 30);

ng = slength(M.gene);

demand_field(M,'V');
nv = slength(M.V);

% convert to Z-scores
Z = nan(ng,nv);
for vi=1:nv
  if length(M.V.val{vi})~=ng, error('V.val must match genelist'); end
  missing = (isnan(M.V.val{vi}) | isinf(M.V.val{vi}));
  mn = mean(M.V.val{vi}(~missing));
  sd = std(M.V.val{vi}(~missing));
  Z(:,vi) = (M.V.val{vi}-mn)./sd;
  Z(missing,vi) = nan;
end

nb = P.bagel_size;
M.gene.bagel = nan(ng,nb);
fprintf('Finding bagels:');
for gi=1:ng, if ~mod(gi,1000), fprintf(' %d/%d',gi,ng); end

  gz = Z(gi,:);
  vidx = find(~isnan(gz));
  if isempty(vidx), continue; end   % no data = no bagel!
  
  others = [1:(gi-1) (gi+1):ng];
  dist2 = bsxfun(@minus,Z(others,vidx),Z(gi,vidx)).^2;
  
  % find nearest neighbors
  method = 2;
  if method == 1
    
    dist2 = nanmean(dist2,2);
    [tmp,ord] = sort(dist2);
    
  elseif method == 2
    % favor genes with more information
    %  each category can contribute at most 10 points; this is downweighted by dist^2
    
    pts = 10*ones(size(dist2));
    pts(dist2>0.001) = 9;
    pts(dist2>0.003) = 8;
    pts(dist2>0.01) = 7;
    pts(dist2>0.03) = 6;
    pts(dist2>0.1) = 5;
    pts(dist2>0.3) = 4;
    pts(dist2>1) = 3;
    pts(dist2>2) = 2;
    pts(dist2>3) = 1;
    pts(dist2>4) = 0;
    pts(isnan(dist2)) = 0;
    
    pts = sum(pts,2);
    [tmp,ord] = sort(pts,'descend');
    
  end
  
  % find nearest neighbors
  bagel = others(ord(1:nb));
  M.gene.bagel(gi,:) = bagel;
end, fprintf('\n');


