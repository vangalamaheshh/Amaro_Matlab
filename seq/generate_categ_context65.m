function run(cset)

build = 'hg18'

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);

if c<1 || c>24, error('c must be 1-24'); end

categ_list.num = generate_categ_context65_names();

d = upper(genome_region(c,1,inf,build));
len = length(d);

% assign "type" category to each position

cen = zeros(1,len);
cen(d=='A') = 1;
cen(d=='C') = 2;
cen(d=='G') = 3;
cen(d=='T') = 4;

clear d;

% assign categories

categ = 16*(cen(2:end-1)-1) + 4*(cen(1:end-2)-1) + cen(3:end);
categ(cen(2:end-1)==0 | cen(1:end-2)==0 | cen(3:end)==0) = 65;
categ = [65 categ 65]';

% SAVE
dir = '/xchip/tcga_scratch/lawrence/db/context65';
if ~exist(dir,'dir'), mkdir(dir); end
if c==1, save_struct(categ_list,[dir '/categs.txt']); end
fname = [dir '/chr' num2str(c) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end % next chr

