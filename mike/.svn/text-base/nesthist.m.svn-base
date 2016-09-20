function nesthist(varargin)

vals = varargin(1:end-1);
nvals = length(vals);
bins = varargin{end};

colors = [0 0 1;0.2 0.7 0.2;1 0 0];
if nvals>3, error('need to define more colors'); end

mn = inf; mx=-inf;
for i=1:nvals
  mn = min(union(mn,vals{i}));
  mx = max(union(mx,vals{i}));
end

if numel(bins)==1
  bins = (mn:((mx-mn)/bins):mx);
end

clf
hold on
for i=1:nvals
  h = histc(vals{i},bins);
  bar(bins,h,1,'facecolor',colors(i,:));
end
xlim([mn mx]);
hold off





