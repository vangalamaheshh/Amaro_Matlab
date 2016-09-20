function out = collapse_strands(in)

c128 = load_struct('/xchip/tcga_scratch/lawrence/db/allcateg/categs128.txt');
c64 = reorder_struct(c128,1:64);

tmp64 = parse(c64.name,'^.*:(.) in (.)_(.)$',{'mid','left','right'});

rc('ACGT')='TGCA';

out = nan(64,size(in,2));
for i=1:64
  j = find(strcmp(tmp64.mid,rc(tmp64.mid{i})) & strcmp(tmp64.left,rc(tmp64.right{i})) &...
            strcmp(tmp64.right,rc(tmp64.left{i})));
  fprintf('%d -> %d\n',i,j);
end
