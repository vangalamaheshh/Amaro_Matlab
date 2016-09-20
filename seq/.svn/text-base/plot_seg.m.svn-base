function Z = plot_seg(S,params)
% plot_seg(S,params)
%
% S is a struct with the following fields:
%     id, chr, start, end, nprobes, segmean
%
% plots it as a vertical red/white/blue plot
%
% S can be a seg array, in which case plots are displayed horizontally in order
%
% Mike Lawrence

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'log',10);
params = impose_default_value(params,'sort_order','none');
params = impose_default_value(params,'xlabels',[]);
params = impose_default_value(params,'rotate_xlabels',true);
params = impose_default_value(params,'max_chr',24);

if isempty(S), error('S empty'); end
if ~iscell(S), S = {S}; end
if any(cellfun('isempty',S)), error('S empty'); end

for i=1:length(S)
  S{i} = require_fields_with_convert(S{i},{'id','chr','start','end','nprobes','segmean'},...
     {{'individual','sample','ID','Sample'},{'chromosome','chrom', 'Chromosome',},{'locstart', 'Start'},{'locend', 'End'},{'notused','nummark', 'Num_Probes'}, ...
      {'log2copyratio','segmean', 'Segment_Mean'}});
  if ~isnumeric(S{i}.chr), S{i}.chr = convert_chr(S{i}.chr); end
  S{i} = make_numeric(S{i},{'start','end','nprobes','segmean'});
end

% load region file

region_list = '/xchip/cga1/lawrence/db/chunks1e6.txt';
R = load_region_file(region_list);
R = reorder_struct(R,R.chr<=params.max_chr);
R2 = keep_fields(R,{'name','chr','start','end'});
R.len = R.end-R.start+1;
% load arm boundaries
R.arm = cell(slength(R),1);
cen = load_cen;
for c=1:params.max_chr
  pq = find(R.chr==c); R.arm(pq) = repmat({[num2str(c) 'p']},length(pq),1);  q = pq(R.start(pq)>cen(c,2)); R.arm(q) = repmat({[num2str(c) 'q']},length(q),1);
end

% load tumor seg file(s) and map to regions

Z = nan(slength(R),length(S));
for chr=1:params.max_chr
  ridx = find(R.chr==chr);
  for s=1:length(S)
    if length(params.log)==1, base = params.log; else base = params.log(s); end
    sidx = find(S{s}.chr==chr);
    for i=1:length(sidx)
      idx = find(R.start(ridx) >= S{s}.start(sidx(i)) & R.end(ridx) <= S{s}.end(sidx(i)));
      Z(ridx(idx),s) = base.^S{s}.segmean(sidx(i));
end,end,end

I = [Z];
I(isnan(I)) = 1;

% sort?

if strcmpi(params.sort_order,'std')
  [tmp ord] = sort(std(I,1));
elseif strcmpi(params.sort_order,'none')
  ord = (1:size(I,2))';
else
  error('unknown params.sort_order = %s',params.sort_order);
end

% figure
clf
imagesc(I(:,ord),[0 2]);
colormap(bwr);colorbar
ylabels_by_group(R.chr);
if isempty(params.xlabels)
  set(gca,'xtick',[]);
else
  set(gca,'xtick',1:length(params.xlabels),'xticklabel',params.xlabels(ord));
  if params.rotate_xlabels
    xticklabel_rotate(1:length(params.xlabels),90,params.xlabels(ord));
  end
end

set(gca,'tickdir','out');
set(gca,'position',[0.08 0.25 0.8 0.7])

if nargout==0
  clear Z
end
