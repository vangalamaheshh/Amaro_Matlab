function [Lt Ln R T N Z]  = load_and_process_CNQC_data(tumor_rcl,normal_rcl,...
  tumor_lanelist,normal_lanelist,lane_blacklist,region_list,normals_db,tumor_seg,normal_seg,params)
% load all data for CopynumberQC analysis
%
% used by fh_CopyNumberQCReport
%
% Mike Lawrence 2009-10-28

if ~exist('tumor_seg','var'), tumor_seg = []; end
if ~exist('normal_seg','var'), normal_seg = []; end

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'min_nprobes',50);
params = impose_default_value(params,'region_median_cutoff',0.2);
params = impose_default_value(params,'choose_random_normals',false);   % (for ruling out weird effects from tumors in normal panel)
params = impose_default_value(params,'ignore_XY',true);
params = impose_default_value(params,'min_lane_nreads',200000);
params = impose_default_value(params,'min_lane_median_nreads',100);
params = impose_default_value(params,'identity_cutoff_euclidean',1e-5);  % for 2000 chunks--values should really be divided by sqrt(n)
params = impose_default_value(params,'collapse_method','none');  % 'gene','segment','10mb', or 'none'
params = impose_default_value(params,'log',2);

% load lanelists
Lt = load_struct(tumor_lanelist);
Ln = load_struct(normal_lanelist);

% identify blacklisted lanes
if strcmpi(lane_blacklist,'none')
  BL = [];
else
  BL = load_lines(lane_blacklist);
end
Lt.is_blacklisted = is_blacklisted(Lt.PU,BL);
Ln.is_blacklisted = is_blacklisted(Ln.PU,BL);

% load region file
R = load_region_file(region_list);
R2 = keep_fields(R,{'name','chr','start','end'});
R.len = R.end-R.start+1;

% load seg file(s) (ground truth), if available
S = cell(2,1);
if ~isempty(tumor_seg)
  if ~exist(tumor_seg,'file')
    fprintf('Not found: %s\n',tumor_seg);
  else
    S{1} = load_cnseg_file(tumor_seg);
    fprintf('Tumor segfile loaded: %s\n',tumor_seg);
    S{1} = reorder_struct(S{1},S{1}.nprobes>=params.min_nprobes);
  end
end
if ~isempty(normal_seg)
  if ~exist(normal_seg,'file')
    fprintf('Not found: %s\n',normal_seg);
  else
    S{2} = load_cnseg_file(normal_seg);
    fprintf('Normal segfile loaded: %s\n',normal_seg);
    S{2} = reorder_struct(S{2},S{2}.nprobes>=params.min_nprobes);
  end
end

% load RCLs (tumor, normal, and panel-of-normals)

fprintf('Loading tumor RCL\n');
[P T] = load_RCL_file(tumor_rcl);
if ~struct_equals(P,R2), error('Tumor RCL doesn''t match region list'); end
if size(T,2) ~= slength(Lt), error('Tumor lanecounts don''t match'); end

fprintf('Loading normal RCL\n');
[P N] = load_RCL_file(normal_rcl);
if ~struct_equals(P,R2), error('Tumor RCL doesn''t match region list'); end
if size(N,2) ~= slength(Ln), error('Normal lanecounts don''t match'); end

fprintf('Loading panel of normals: ');
pnlist = setdiff(load_lines(normals_db),{''});

Q = cell(length(pnlist),1);
for i=1:length(pnlist),if ~mod(i,100), fprintf('%d/%d ',i,length(pnlist)); end
  try
    [P Q{i}] = load_RCL_file(pnlist{i});
    if ~struct_equals(P,R2)
      fprintf('Normal panel RCL #%d doesn''t match region list\n',i);
      Q{i}=[];
      continue
    end
  catch me
    fprintf('Error loading normal panel RCL #%d: %s\n',i,me.message);
    Q{i}=[];
    continue
  end
end,fprintf('\n');
Q = cat(2,Q{:});

if size(Q,2)<5, error('Panel of normals has %d lanes: need at least 5.',size(Q,2)); end

% ignore X and Y chromosomes
if params.ignore_XY
  xy = find(R.chr>22);
  T(xy,:)=[]; N(xy,:)=[]; Q(xy,:)=[];
  R = reorder_struct(R,setdiff(1:slength(R),xy));
end

% collapse data (if requested)
switch params.collapse_method

 case 'none'
  % (use this option if capture data is already gathered according to chunk1e7)

 case 'gene'
  % collapse RCL data to "genes"
  %   for capture data: collapses regions to genes
  %   for WGS data: collapses chunks to themselves (no processing necessary)
  [g gi gj] = unique(R.name);
  if length(g) < slength(R)   % if any collapsing needs to happen
    fprintf('Collapsing all data to genes... ');
    for i=1:length(g),if ~mod(i,1000), fprintf('%d/%d ',i,length(g)); end
      idx = find(gj==i);
      T(gi(i),:) = sum(T(idx,:),1);
      N(gi(i),:) = sum(N(idx,:),1);
      Q(gi(i),:) = sum(Q(idx,:),1);
      R.start(gi(i)) = min(R.start(idx));
      R.end(gi(i)) = max(R.end(idx));
      R.len(gi(i)) = sum(R.len(idx));
    end, fprintf('\n');
    [R ord] = sort_struct(reorder_struct(R,gi),{'chr','start'});   % re-sort to genome order
    T = T(gi(ord),:); N = N(gi(ord),:); Q = Q(gi(ord),:);
  end

 case 'segment'
  % collapse RCL data to tumor segments (from array)
  if isempty(S{1}), error('No tumor segment data available'); end
  fprintf('Collapsing all data to %d tumor segments\n',slength(S{1}));
  To=T; No=N; Qo=Q; Ro=R;
  R = keep_fields(S{1},{'chr','start','end'});
  T = nan(slength(R),size(To,2));
  N = nan(slength(R),size(No,2));
  Q = nan(slength(R),size(Qo,2));
  for i=1:slength(R)
    idx = find(Ro.chr==R.chr(i) & Ro.start>=R.start(i) & Ro.end<=R.end(i));
    T(i,:) = sum(To(idx,:),1);
    N(i,:) = sum(No(idx,:),1);
    Q(i,:) = sum(Qo(idx,:),1);
  end

 case '10mb'
  % collapse RCL data to 10Mb chunks
  fprintf('Collapsing all data to 10Mb chunks\n');
  To=T; No=N; Qo=Q; Ro=R;
  R = load_region_file('/xchip/tcga_scratch/lawrence/db/chunks1e7.txt');
  R2 = keep_fields(R,{'name','chr','start','end'});
  R.len = R.end-R.start+1;
  T = nan(slength(R),size(To,2));
  N = nan(slength(R),size(No,2));
  Q = nan(slength(R),size(Qo,2));
  for i=1:slength(R)
    idx = find(Ro.chr==R.chr(i) & Ro.start>=R.start(i) & Ro.end<=R.end(i));
    T(i,:) = sum(To(idx,:),1);
    N(i,:) = sum(No(idx,:),1);
    Q(i,:) = sum(Qo(idx,:),1);
  end

 otherwise
  error('Unknown params.collapse_method');
end

% exclude lanes with too few counts or very low median counts
qcov = sum(Q,1)'; Lt.cov = sum(T,1)'; Ln.cov = sum(N,1)';
mq = median(Q,1)';
totok = (qcov>params.min_lane_nreads);
medok = (mq>params.min_lane_median_nreads);
%qok = totok & medok;
qok = totok;
Ln.enoughreads = (Ln.cov>params.min_lane_nreads);
Lt.enoughreads = (Lt.cov>params.min_lane_nreads);
Q = Q(:,qok);

if ~any(Ln.enoughreads)
  fprintf('WARNING: No normal lanes have enough reads\n');
  Ln.use = true(slength(Ln),1);
else
  Ln.use = Ln.enoughreads;
end  

if ~any(Lt.enoughreads)
  fprintf('WARNING: No tumor lanes have enough reads\n');
  Lt.use = true(slength(Lt),1);
else
  Lt.use = Lt.enoughreads;
end

% normalize all lanes to total lane counts
Q = bsxfun(@rdivide,Q,qcov(qok)');
T = bsxfun(@rdivide,T,Lt.cov');
N = bsxfun(@rdivide,N,Ln.cov');

% for each test-normal, find 5-nearest-normals in the panel
% (any panel-normal that is an "exact match" to any normal lane, is excluded from the panel entirely)
fprintf('Finding nearest normals\n');
d = dist(N',Q','euclidean');
notexactmatches = find(all(d(Ln.use,:)>params.identity_cutoff_euclidean,1));
Q = Q(:,notexactmatches);
d = d(:,notexactmatches);
if size(Q,2)<5, error('Reduced panel of normals has %d lanes: need at least 5.',size(Q,2)); end
[tmp ord] = sort(d,2);
nn = ord(:,1:5);

if params.choose_random_normals
  % (for ruling out weird effects from tumors or tumor-in-normals in the panel of normals)
  fprintf('\n***********************\nCHOOSING RANDOM NORMALS\n**********************\n');
  nn = ceil(size(Q,2)*rand(size(nn)))
  %nn = ceil(5*rand(size(nn)))
end

% divide each normal lane by the median of its 5-nearest-normals
for i=1:size(N,2)
  N(:,i) = N(:,i) ./ median(Q(:,nn(i,:)),2);
end

% reduce panel of normals to the union of all 5-nearest-normals for the test normals
unn = unique(nn(Ln.use,:));
Q = Q(:,unn);
d = d(:,unn);

% for each test-tumor, find 5-nearest-normals in the reduced panel
d = dist(T',Q','euclidean');
[tmp ord] = sort(d,2);
nn = ord(:,1:5);

% divide each tumor lanes by the median of its 5-nearest-normals
for i=1:size(T,2)
  T(:,i) = T(:,i) ./ median(Q(:,nn(i,:)),2);
end

% get rid of regions with very low coverage across the panel of normals
m = median(Q,2);
keep = find(m > median(m)*params.region_median_cutoff);
R = reorder_struct(R,keep);
Q = Q(keep,:); T = T(keep,:); N = N(keep,:);

% map seg file(s), if available, to final region list
Z = nan(slength(R),2);
for s=1:2
  if ~isempty(S{s})
    for r=1:slength(R)
      sidx = find(S{s}.chr==R.chr(r) & S{s}.start<=R.end(r) & S{s}.end>=R.start(r));
      % weight by amount of overlap
      st = max(R.start(r),S{s}.start(sidx));
      en = min(R.end(r),S{s}.end(sidx));
      Z(r,s) = weighted_mean(params.log.^S{s}.segmean(sidx),en-st+1);
    end
  end
end

% compute noise
Lt.noise = median(abs(diff([T.^3])))';
Ln.noise = median(abs(diff([N.^3])))';

% compute normalness >0.9 = normal, 0.8-0.9=uncertain, <0.8=tumor
Lt.normalness = 1-std(log2(T))';
Ln.normalness = 1-std(log2(N))';

% measure correlations
%   (1) to tumor seg (if available)
%   (2) to median of tumor lanes

have_tumor_seg = ~all(isnan(Z(:,1)));
if have_tumor_seg
  tumorseg = Z(:,1);
  good = find(~isnan(tumorseg));
  Lt.tumorseg_corr = 1 - dist(T(good,:)',tumorseg(good)','correlation');
  Ln.tumorseg_corr = 1 - dist(N(good,:)',tumorseg(good)','correlation');
else
  Lt.tumorseg_corr = nan(size(T,2),1);
  Ln.tumorseg_corr = nan(size(N,2),1);
end

tumormed = median(T(:,Lt.use),2);
good = find(~isnan(tumormed));
Lt.tumormed_corr = 1 - dist(T(good,:)',tumormed(good)','correlation');
Ln.tumormed_corr = 1 - dist(N(good,:)',tumormed(good)','correlation');

% ASSIGN JUDGEMENTS

% normal lanes:
%   PROBABLY_OK              unless:
%   SUSPICIOUSLY_ABNORMAL    if tumor(seg|med)_corr<0.6 and normalness<0.85
%   MIXUP_UNKNOWN_TUMOR      if [tumor(seg|med)_corr<0.6 and normalness<0.75] or normalness==nan      -->*mixup*
%   MIXUP_MATCHED_TUMOR      if tumor(seg|med)_corr>=0.6 & normalness<mxntl+0.1    -->*mixup*
%   CONTAM_TUMOR_IN_NORMAL   if tumor(seg|med)_corr>=0.6 & normalness>mxntl+0.1    -->*mixup*
%                            (mxntl = max normalness of tumor lanes)
%   LOW_COUNTS               if counts<threshold

% tumor lanes:
%   PROBABLY_OK              unless:
%   IDENTITY_CONFIRMED       if tumorseg_corr>=0.7
%   MIXUP_WRONG_TUMOR        if tumor(seg|med)_corr<0.3                        -->*mixup*
%   MIXUP_NORMAL             if normalness>=0.9                                -->*mixup*
%   LOW_COUNTS               if counts<threshold

% NORMAL judgements

Ln.judgement = repmat({'PROBABLY_OK'},slength(Ln),1);
idx = find(Ln.normalness<0.85); Ln.judgement(idx) = repmat({'SUSPICIOUSLY_ABNORMAL'},length(idx),1);
idx = find(Ln.normalness<0.75 | isnan(Ln.normalness)); Ln.judgement(idx) = repmat({'MIXUP_UNKNOWN_TUMOR'},length(idx),1);
mnntl = min(Lt.normalness);
mxntl = max(Lt.normalness);
if have_tumor_seg, tumor_corr_N = Ln.tumorseg_corr; else tumor_corr_N = Ln.tumormed_corr; end
if mnntl<=0.8   % check if the normal has the matched tumor (skip if tumor is extremely CN quiet)
  idx = find(Ln.normalness<mxntl+0.1 & tumor_corr_N>=0.6); Ln.judgement(idx) = repmat({'MIXUP_MATCHED_TUMOR'},length(idx),1);
  idx = find(Ln.normalness>mxntl+0.1 & tumor_corr_N>=0.6); Ln.judgement(idx) = repmat({'CONTAM_TUMOR_IN_NORMAL'},length(idx),1);
end
idx = find(~Ln.enoughreads); Ln.judgement(idx) = repmat({'LOW_COUNTS'},length(idx),1);

% TUMOR judgements

Lt.judgement = repmat({'PROBABLY_OK'},slength(Lt),1);
idx = find(Lt.tumorseg_corr>=0.7); Lt.judgement(idx) = repmat({'IDENTITY_CONFIRMED'},length(idx),1);
if have_tumor_seg, tumor_corr_T = Lt.tumorseg_corr; else tumor_corr_T = Lt.tumormed_corr; end
idx = find(tumor_corr_T<0.3); Lt.judgement(idx) = repmat({'MIXUP_WRONG_TUMOR'},length(idx),1);
idx = find(Lt.normalness>=0.9); Lt.judgement(idx) = repmat({'SUSPICIOUSLY_NORMAL'},length(idx),1);

% mark noisy lanes (supersedes "CONTAM/MIXUP" judgement)
idx = find(Lt.noise>=0.4); Lt.judgement(idx) = repmat({'EXTREMELY_NOISY'},length(idx),1);
idx = find(Ln.noise>=0.4); Ln.judgement(idx) = repmat({'EXTREMELY_NOISY'},length(idx),1);

% mark low-coutns lanes (supersedes everything)
idx = find(~Lt.enoughreads); Lt.judgement(idx) = repmat({'LOW_COUNTS'},length(idx),1);

% flag mixups to report to Firehose
Lt.is_mixup = (contains(Lt.judgement,'MIXUP') | contains(Lt.judgement,'CONTAM'));
Ln.is_mixup = (contains(Ln.judgement,'MIXUP') | contains(Ln.judgement,'CONTAM'));

