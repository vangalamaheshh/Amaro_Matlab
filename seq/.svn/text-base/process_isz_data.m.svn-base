function X = process_isz_data(X,params)
% Mike Lawrence 2009-10-15

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'low_counts_threshold',100000);
params = impose_default_value(params,'width_pctdev_threshold',5);
params = impose_default_value(params,'adjmean_pctdev_threshold',1.2);
params = impose_default_value(params,'pardoned_early_libraries',{'Solexa-9039','Solexa-9040','Solexa-9041'});
                        %  GBM-0188 had some early libraries with highly variable peak width

X.mx = size(X.dat,1)-1;
X.nlanes = size(X.dat,2);
X.mode=nan(X.nlanes,1);X.adjmean=X.mode;X.width=X.mode;
fprintf('Lanes: ');for i=1:X.nlanes, fprintf('%d/%d ',i,X.nlanes);
  X.mode(i) = weighted_mode([0:X.mx]',X.dat(:,i));
  radius = round(2*sqrt(X.mode(i)));
  range = [max(0,X.mode(i)-radius):X.mode(i)+radius]';
  X.adjmean(i) = weighted_mean(range,X.dat(range+1,i));
  if X.mode(i)>0
    q = interp(X.dat(:,i),100);
    [qmax qmode] = max(q);
    qhalfmax = qmax/2;
    left = find(flipud(q(1:qmode))<=qhalfmax,1);
    if isempty(left), left = qmode; end
    right = find(q(qmode:end)<=qhalfmax,1);
    X.width(i) = (right+left)/100;
  end
end,fprintf('\n');
X.cov = sum(X.dat,1)';
X.good = find(X.cov>=params.low_counts_threshold);

% normalize each lane
X.datnorm = X.dat ./ repmat(X.cov',X.mx+1,1);

% compute per-library stats and find outliers

X.lib = [];
[X.lib.name li lj] = unique(X.lanes.LB);
X.nlib = slength(X.lib);
X.lib.median_adjmean = nan(X.nlib,1);
X.lib.median_width = nan(X.nlib,1);
X.adjmean_pctdev = nan(X.nlanes,1);
X.width_pctdev = nan(X.nlanes,1);
X.fewlanes = [];
X.pardoned = [];
for i=1:X.nlib
  idx = find(lj==i);
  idx = intersect(idx,X.good);      % don't look at low-counts lanes
  if length(idx)<3
    X.fewlanes = [X.fewlanes; as_column(idx)];
    continue;
  end         % don't look at libraries with fewer than 3 good lanes
  if ismember(X.lib.name{i},params.pardoned_early_libraries)
    X.pardoned = [X.pardoned; as_column(idx)];
  end
  x = X.adjmean(idx);
  y = X.width(idx);
  mx = median(x);
  my = median(y);
  X.lib.median_adjmean(i) = mx;
  X.lib.median_width(i) = my;
  X.adjmean_pctdev(idx) = 100*(x-mx)/mx;
  X.width_pctdev(idx) = 100*(y-my)/my;
end
X.outliers = find(abs(X.adjmean_pctdev)>params.adjmean_pctdev_threshold | ...
                  abs(X.width_pctdev)>params.width_pctdev_threshold);
X.outliers = setdiff(X.outliers,X.pardoned);

X.blacklisted = find(X.lanes.is_blacklisted);
X.blacklisted_outliers = intersect(X.blacklisted,X.outliers);
X.outliers = setdiff(X.outliers,X.blacklisted_outliers);

