function run(cset)
% each position tells how far away is the nearest C in a "AID hotspot" WRCY sequence (on either strand)
% (from zero, up up to a maximum of 1000 bases); 1001 means the site is an N

build = 'hg18'

dir = '/xchip/cga1/lawrence/db/wrcy';
if ~exist(dir,'dir'), mkdir(dir); end

maxdist = 1000;

if ~exist('cset','var'), cset=1:24; end
for csetidx=1:length(cset)
chr = cset(csetidx);
fprintf('\nchr%d\n',chr);
if chr<1 || chr>24, error('chr must be 1-24'); end

if chr==1
  categ_list = [];
  categ_list.num = (0:maxdist)';
  categ_list.name = num2cellstr(categ_list.num);
  categ_list.num(end+1) = maxdist+1;
  categ_list.name{end+1} = 'base_is_an_N';
  save_struct(categ_list,[dir '/categs.txt']);
end

d = upper(genome_region(chr,1,inf,build));
len = length(d);

% find hotspots
hot = false(len,1);
for i=3:len-2
  if d(i)=='C'
    if d(i-1)=='G' || d(i-1)=='A'
      if d(i-2)=='A' || d(i-2)=='T'
        if d(i+1)=='C' || d(i+1)=='T'
          hot(i)=true;
    end,end,end
  elseif d(i)=='G'
    if d(i-1)=='A' || d(i-1)=='G'
      if d(i+1)=='C' || d(i+1)=='T'
        if d(i+2)=='A' || d(i+2)=='T'
          hot(i)=true;
end,end,end,end,end
hidx = find(hot);

% find distances from nearest hotspot on the left
left = nan(len,1);
dist = maxdist;
for i=1:len
  if hot(i), dist=0; else if dist<maxdist, dist=dist+1; end, end
  left(i)=dist;
end
% find distances from nearest hotspot on the right
right = nan(len,1);
dist = maxdist;
for i=len:-1:1
  if hot(i), dist=0; else if dist<maxdist, dist=dist+1; end, end
  right(i)=dist;
end
% find distances from nearest hotspot on either side
categ = min(left,right);
categ(d=='N') = maxdist+1;

% SAVE
fname = [dir '/chr' num2str(chr) '.mat'];
save(fname,'categ');
f = fopen(regexprep(fname,'\.mat$','\.txt'),'wt');
fprintf(f,'%d\n',categ);
fclose(f);

end % next chr

