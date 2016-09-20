function X = dRanger_get_counts(X,sample,P)
%
% calls extract_from_cbb and extract_from_mqz for each region in the list of rearrangements
%
% Mike Lawrence 2009-06-23

if ~exist('sample','var'), error('<sample> is required'); end

if ~exist('P','var'), P=[]; end
if ~isempty(P) && ~isstruct(P), error('Second parameter should be a P struct'); end

P=impose_default_value(P,'getcounts_tempdir','dRgg');

try

require_fields(X,{'chr1','min1','max1','chr2','min2','max2'});

direc = ['/xchip/tcga_scratch/lawrence/' sample];
tempdir = [direc '/' P.getcounts_tempdir];
if ~exist(tempdir), mkdir(tempdir); end

% output target list
targfile = [tempdir '/targs.txt'];
for i=1:slength(X);tmp{i,1}=num2str(X.num(i));end
T.name = [regexprep(tmp,'(.*)','ra$1end1');regexprep(tmp,'(.*)','ra$1end2')];
T.chr = [X.chr1; X.chr2];
T.start = [X.min1; X.min2];
T.end = [X.max1; X.max2];
save_struct(T,targfile,'no_headers');

% call extract_from_cbb

outfile = [tempdir '/targ_cbb.txt'];
extract_from_cbb(sample,targfile,outfile);

keyboard

catch me; excuse(me); end

