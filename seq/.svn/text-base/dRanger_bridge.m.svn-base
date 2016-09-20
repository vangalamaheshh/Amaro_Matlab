function dRanger_bridge(sample,rnum,P)
% dRanger_bridge(sample,rnum,P)
%
% sample
%
% rnum = which rearrangement (-->"num" column of dRanger_results.txt)
%        [in future version, vector rnum will be allowed]
%
% Mike Lawrence 2009-09-30

if ~exist('sample','var'), error('<sample> is required'); end
if ~exist('rnum','var'), error('<rnum> is required'); end
if ~isnumeric(rnum) || numel(rnum)~=1, error('<rnum> must be a single number'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'results_name','dRanger_results');
P=impose_default_value(P,'radius_to_look',100);

fprintf('\ndRanger_bridge\n\tsample = %s\n\trnum = %d\n\n',sample,rnum);

% directories and filenames

short_sample = sample_to_short_sample(sample);
direc = ['/xchip/tcga_scratch/lawrence/' sample];
fname = [direc '/' P.results_name '.txt'];

% load dRanger results

fprintf('Loading dRanger results\n');
if ~exist(fname,'file'), error('Can''t find file %s',fname);end
X = load_struct(fname);
X = reorder_struct(X,strcmp(X.num,num2str(rnum)));
if length(X)~=1, error('Problem finding rnum %d\n',rnum); end
X = make_numeric(X,{'chr1','str1','pos1','min1','max1','chr2','str2','pos2','min2','max2'});
E = cell(2,1);
E{1}.chr = X.chr1; E{1}.str = X.str1; E{1}.pos = X.pos1;
E{2}.chr = X.chr2; E{2}.str = X.str2; E{2}.pos = X.pos2;

% retrieve reads from EOOM

for e=1:2
  E{e}.left = E{e}.pos - P.radius_to_look;
  E{e}.right = E{e}.pos + P.radius_to_look;
  [E{e}.R E{e}.B E{e}.S E{e}.ri E{e}.bi E{e}.si] = ...
    pull_from_boom([direc '/tumor.boom'],E{e}.chr,E{e}.left,E{e}.right);
  E{e}.nr = size(E{e}.R,1);
  E{e}.nb = size(E{e}.B,1);
  E{e}.ns = size(E{e}.S,1);
end

% for each read, look for a clear split between ref and nonref
% (for str=0, want ref->nonref; for str=1, want nonref->ref)

keyboard

for e=1:2
  xleft = E{e}.left - 200;
  xright = E{e}.right + 200;
  offset = xleft-1;
  x = zeros(xright-offset,1);
  bidx1 = 1;
  for r=1:E{e}.nr
    st = E{e}.R(r,4);
    en = E{e}.R(r,5);
    bidx2 = bidx1 + en -st; 
    base = E{e}.B(bidx1:bidx2,1);
    isnonref = (base>64);
    qual = E{e}.B(bidx1:bidx2,2);
    nonrefqual = isnonref .* qual;
    d = zeros(length(nonrefqual),1);
    for i=2:length(nonrefqual)-1
      d(i) = mean(nonrefqual(i+1:end))-mean(nonrefqual(1:i));
    end
    x(st-offset:en-offset) = x(st-offset:en-offset) + d;
    bidx1 = bidx2 + 1;
  end
  if E{e}.str==0, [tmp idx] = max(x); else [tmp idx] = min(x); end
  pos = idx+offset;
  fprintf('end%d: breakpoint = %d\n',e,pos);

  %plot(x)
end
  

















% throw out all-reference reads

min_nonrefqual = 10;
min_nonrefbases = 2;
for e=1:2
  offset = E{e}.left-1;
  bidx = find(E{e}.E(:,1)>=65 & E{e}.E(:,2)>=min_nonrefqual);
%  h = histc(E{e}.E(bidx,3)
  [r ri rj] = unique(E{e}.E(bidx,3));
  
  keyboard
end



% look for non-ref consensus

minqual = 10;
for e=1:2
  offset = E{e}.left-1;
  E{e}.En = E{e}.E(E{e}.E(:,1)>=65 & E{e}.E(:,2)>=minqual,:);
  E{e}.qn = zeros(4,E{e}.right-offset);
  for pos=E{e}.left:E{e}.right
    idx = find(E{e}.En(:,4)==pos);
    for i=1:4
      E{e}.qn(i,pos-offset) = sum(E{e}.En(idx(E{e}.En(idx,1)==i+64),2));
    end
  end
end

keyboard

% slide the regions against each other, compute a score for each alignment

for r=-P.radius_to_look:P.radius_to_look
end  




plot(smooth(max(E{2}.qn),10))
