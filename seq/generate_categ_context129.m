function run(cset)
% generates context file for expanded context list
%
% output: list of contexts across all genome positions
%
% contexts 1-64 =   trinucleotide identity (in transcribed region, with respect to transcribed strand)
% contexts 65-128 = trinucleotide identity (in untranscribed region, with respect to (+) strand)
% context 129 =     all N's, bases directly adjacent to N's, and the first and last bases in the chromosome.
%
% Mike Lawrence 2009-10-27

build = 'hg18'

R = load_refseq(build);

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);

if c<1 || c>24, error('c must be 1-24'); end

categ_list.num = (1:129)';
x = {'A';'C';'G';'T'};
y = {}; for i=1:length(x), y = [y; regexprep(x,'^(.*)$',[x{i} '_$1'])]; end
z = {}; for i=1:length(x), z = [z; regexprep(y,'^(.*)$',[x{i} ' in $1'])]; end
categ_list.name = [regexprep(z,'^(.*)$',['transcribed $1']);regexprep(z,'^(.*)$',['nontranscribed $1']);'any N'];

d = upper(genome_region(c,1,inf,build));
len = length(d);

fprintf('Annotating...\n');

% find out which positions are transcribed, and from what strand (i.e. in target list)
Tc = reorder_struct(T,T.chr==c & T.start<=len);
Tc.end = min(len,Tc.end);
tx = nan(len,1);
% nan = nontranscribed; 0 = (+)strand transcribed; 1 = (-)strand transcribed
for i=1:slength(Tc), tx(Tc.start(i):Tc.end(i)) = Tc.strand(i); end

% assign "type" category to each position

cen = zeros(1,len);
cen(d=='A') = 1;
cen(d=='C') = 2;
cen(d=='G') = 3;
cen(d=='T') = 4;

clear d;

left = [0 cen(1:end-1)];
right = [cen(2:end) 0];

% assign categories

categ = nan(1,len);

% nontranscribed
idx = find(isnan(tx));
t = 64 + (16*(cen(idx)-1) + 4*(left(idx)-1) + right(idx));
categ(idx) = t;

% (+) transcribed
idx = find(tx==0);
t = 16*(cen(idx)-1) + 4*(left(idx)-1) + right(idx);
categ(idx) = t;

% (-) transcribed
idx = find(tx==1);
rc = [4 3 2 1];
t = 16*(rc(cen(idx))-1) + 4*(rc(right(idx))-1) + rc(left(idx));
categ(idx) = t;

categ(cen==0 | left==0 | right==0) = 129;

% SAVE
dir = '/xchip/tcga_scratch/lawrence/db/context129';
if ~exist(dir,'dir'), mkdir(dir); end
if c==1, save_struct(categ_list,[dir '/categs.txt']); end
fname = [dir '/chr' num2str(c) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end % next chr

