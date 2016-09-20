function run(cset)
% 1-1024 = AAAAA-TTTTT
% 1025 = N or end within 2 bases

build = 'hg18'

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
c = cset(csetidx);
fprintf('\nchr%d\n',c);

if c<1 || c>24, error('c must be 1-24'); end

categ_list=[];
categ_list.num = (1:1024)';
b = 'ACGT';
bi = nan(1024,5);
for i=1:1024
  t=i-1; for j=1:5, bi(i,j) = mod(t,4)+1; t=floor(t/4); end
  categ_list.name{i,1} = [b(bi(i,5)) ' in ' b(bi(i,[1 3])) '_' b(bi(i,[4 2]))];
end
categ_list.num(end+1) = 1025;
categ_list.name{end+1} = 'any N';

% load chromosome reference sequence
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

categ = 256*(cen(3:end-2)-1)+...
         64*(cen(4:end-1)-1)+...
         16*(cen(2:end-3)-1)+...
          4*(cen(5:end-0)-1)+...
             cen(1:end-4); 

categ(cen(1:end-4)==0 | cen(2:end-3)==0 | cen(3:end-2)==0 | cen(4:end-1)==0 | cen(5:end)==0) = 1025;
categ = [1025 1025 categ 1025 1025]';

% SAVE
dir = '/xchip/cga1/lawrence/db/pentamer';
if ~exist(dir,'dir'), mkdir(dir); end
if c==1, save_struct(categ_list,[dir '/categs.txt']); end
fname = [dir '/chr' num2str(c) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end % next chr

