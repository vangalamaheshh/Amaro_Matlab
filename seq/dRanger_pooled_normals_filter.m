function dRanger_pooled_normals_filter(samples,P)
% dRanger_pooled_normals_filter(samples,P)
%
% OBSOLETE:
%    REPLACED BY    dRanger_screen_panel_of_normals
%
% given a set of dRanger results files
%   -- specified by a list of samples and P.results_name
%
% updates the columns "somatic_in_other_samples" and "germline_in_other_samples"
%   -- tells how many of the _other_ results files have an entry for the
%      same* rearrangement with normreads==0 (-->"other_samples_somatic")  == recurrent somatic
%                            or normreads>0  (-->"other_samples_germline") == probable artifact/polymorphism
%
% Note: Each other sample contributes at most 1 to "other_samples_somatic" OR "other_samples_germline", NOT BOTH.
%       In cases where there is a conflict, germline status takes precedence.
%
% also updates the column "num_other_samples" = number of other results files examined.
%
% *same = chr1, chr2, str1, str2 are same; (min1-max1) overlaps; (min2-max2) overlaps
%
% Mike Lawrence 2009-10-13

if ~exist('samples','var'), error('<samples> is required'); end
if ~iscell(samples), error('<samples> should be a list of samples (cell array)'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P = impose_default_value(P,'results_name','dRanger_results');
P = impose_default_value(P,'extra_match_margin',0);

basedir = '/xchip/tcga_scratch/lawrence';
ns = length(samples);

fprintf('Checking for availability of all results files... ');
fname = cell(ns,1);
for s=1:ns
  fname{s} = [basedir '/' samples{s} '/' P.results_name '.txt'];
  if ~exist(fname{s},'file'), error('\nNot found: %s',fname{s}); end
end, fprintf('OK.\n');

fprintf('Loading results for all samples: ');
X = cell(ns,1);
nr = zeros(ns,1);
flds = {'chr1','chr2','str1','str2','min1','min2','max1','max2','normreads','filterB','filterQ','filterW'};
for s=1:ns, fprintf('%d/%d ',s,ns);
  X{s} = load_struct(fname{s});
  require_fields(X{s},flds);
  X{s} = make_numeric(X{s},flds);
  nr(s) = slength(X{s});
end, fprintf('\n');

% build unified list

CHR1=1; STR1=2; MIN1=3; MAX1=4;
CHR2=5; STR2=6; MIN2=7; MAX2=8;
NORMREADS=9; FILTERB=10; FILTERQ=11; FILTERW=12;
SAMPLE=13; ROW=14;
OTHER_SOM=15; OTHER_GERM=16;

x = nan(sum(nr),OTHER_GERM);
pos1 = 1;
for s=1:ns
  pos2 = pos1+nr(s)-1;
  x(pos1:pos2,CHR1) = X{s}.chr1;
  x(pos1:pos2,STR1) = X{s}.str1;
  x(pos1:pos2,MIN1) = X{s}.min1;
  x(pos1:pos2,MAX1) = X{s}.max1;
  x(pos1:pos2,CHR2) = X{s}.chr2;
  x(pos1:pos2,STR2) = X{s}.str2;
  x(pos1:pos2,MIN2) = X{s}.min2;
  x(pos1:pos2,MAX2) = X{s}.max2;
  x(pos1:pos2,NORMREADS) = X{s}.normreads;
  x(pos1:pos2,FILTERB) = X{s}.filterB;
  x(pos1:pos2,FILTERQ) = X{s}.filterQ;
  x(pos1:pos2,FILTERW) = X{s}.filterW;
  x(pos1:pos2,SAMPLE) = s*ones(nr(s),1);
  x(pos1:pos2,ROW) = (1:nr(s))';
  pos1 = pos2+1;
end
x = sortrows(x);
nx = size(x,1);

% traverse list
%   maintain pointers:
%
%   upper = the topmost record whose end1 overlaps record "idx"
%   idx   = the record we're currently adjudicating 
%   lower = the bottommost record whose end1 overlaps record "idx"

x(:,[OTHER_SOM OTHER_GERM]) = nan;
upper = 1; lower = 1;
for idx=1:nx, if ~mod(idx,10000), fprintf('%d/%d ',idx,nx); end
  while x(upper,CHR1)<x(idx,CHR1) || x(upper,STR1)<x(idx,STR1) || x(upper,MAX1)<x(idx,MIN1)-P.extra_match_margin, upper=upper+1; end
  while x(lower,CHR1)<x(idx,CHR1) || x(lower,STR1)<x(idx,STR1) || ...
    (lower<nx && x(lower+1,CHR1)==x(idx,CHR1) && x(lower+1,STR1)==x(idx,STR1) && ...
    x(lower+1,MIN1)-P.extra_match_margin<=x(idx,MAX1)), lower=lower+1; end
  match = upper:lower;
  match = match(x(match,CHR2)==x(idx,CHR2));
  match = match(x(match,STR2)==x(idx,STR2));
  match = match(x(match,SAMPLE)~=x(idx,SAMPLE));
  match = match(x(match,MIN2)-P.extra_match_margin<=x(idx,MAX2));
  match = match(x(match,MAX2)>=x(idx,MIN2)-P.extra_match_margin);
  germ = match(x(match,NORMREADS)>0);
  germ_samps = unique(x(germ,SAMPLE)); 
  som_samps = setdiff(x(match,SAMPLE),germ_samps);
  x(idx,OTHER_SOM) = length(som_samps);
  x(idx,OTHER_GERM) = length(germ_samps);
end, fprintf('\n');








% print out some filtering stats

som = find(x(:,NORMREADS)==0);
som_q = som(x(som,FILTERQ)==0);
som_w = som(x(som,FILTERW)==0);
som_b = som(x(som,FILTERB)==0);
som_qw = intersect(som_q,som_w);
som_bw = intersect(som_b,som_w);
som_bq = intersect(som_b,som_q);
som_bqw = intersect(som_b,som_qw);

germ = find(x(:,NORMREADS)>0);
germ_q = germ(x(germ,FILTERQ)==0);
germ_w = germ(x(germ,FILTERW)==0);
germ_b = germ(x(germ,FILTERB)==0);
germ_qw = intersect(germ_q,germ_w);
germ_bw = intersect(germ_b,germ_w);
germ_bq = intersect(germ_b,germ_q);
germ_bqw = intersect(germ_b,germ_qw);

o00 = find(x(:,OTHER_SOM)==0 & x(:,OTHER_GERM)==0);
o10 = find(x(:,OTHER_SOM)>=1 & x(:,OTHER_GERM)==0);
o11 = find(x(:,OTHER_SOM)>=1 & x(:,OTHER_GERM)>=1);
o01 = find(x(:,OTHER_SOM)==0 & x(:,OTHER_GERM)>=1);

m_som = [...
length(som) length(intersect(som,o00)) length(intersect(som,o10)) length(intersect(som,o11)) length(intersect(som,o01));...
length(som_b) length(intersect(som_b,o00)) length(intersect(som_b,o10)) length(intersect(som_b,o11)) length(intersect(som_b,o01));...
length(som_q) length(intersect(som_q,o00)) length(intersect(som_q,o10)) length(intersect(som_q,o11)) length(intersect(som_q,o01));...
length(som_w) length(intersect(som_w,o00)) length(intersect(som_w,o10)) length(intersect(som_w,o11)) length(intersect(som_w,o01));...
length(som_bq) length(intersect(som_bq,o00)) length(intersect(som_bq,o10)) length(intersect(som_bq,o11)) length(intersect(som_bq,o01));...
length(som_bw) length(intersect(som_bw,o00)) length(intersect(som_bw,o10)) length(intersect(som_bw,o11)) length(intersect(som_bw,o01));...
length(som_qw) length(intersect(som_qw,o00)) length(intersect(som_qw,o10)) length(intersect(som_qw,o11)) length(intersect(som_qw,o01));...
length(som_bqw) length(intersect(som_bqw,o00)) length(intersect(som_bqw,o10)) length(intersect(som_bqw,o11)) length(intersect(som_bqw,o01));...
]

 
m_germ = [...
length(germ) length(intersect(germ,o00)) length(intersect(germ,o10)) length(intersect(germ,o11)) length(intersect(germ,o01));...
length(germ_b) length(intersect(germ_b,o00)) length(intersect(germ_b,o10)) length(intersect(germ_b,o11)) length(intersect(germ_b,o01));...
length(germ_q) length(intersect(germ_q,o00)) length(intersect(germ_q,o10)) length(intersect(germ_q,o11)) length(intersect(germ_q,o01));...
length(germ_w) length(intersect(germ_w,o00)) length(intersect(germ_w,o10)) length(intersect(germ_w,o11)) length(intersect(germ_w,o01));...
length(germ_bq) length(intersect(germ_bq,o00)) length(intersect(germ_bq,o10)) length(intersect(germ_bq,o11)) length(intersect(germ_bq,o01));...
length(germ_bw) length(intersect(germ_bw,o00)) length(intersect(germ_bw,o10)) length(intersect(germ_bw,o11)) length(intersect(germ_bw,o01));...
length(germ_qw) length(intersect(germ_qw,o00)) length(intersect(germ_qw,o10)) length(intersect(germ_qw,o11)) length(intersect(germ_qw,o01));...
length(germ_bqw) length(intersect(germ_bqw,o00)) length(intersect(germ_bqw,o10)) length(intersect(germ_bqw,o11)) length(intersect(germ_bqw,o01));...
]
