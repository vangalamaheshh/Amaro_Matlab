%%%%   regulatory regions
%%%%
%%%%   (defined by 3-mammal alignment from Mike Chapman)

load('/xchip/cga1/lawrence/mm/analysis/20100318_cons/regulatory_regions.mat','R');

chrlen = load_chrlen;
outd = '/xchip/cga1/lawrence/db/regulatory';
mkdir(outd);

for c=1:24,disp(c);
  len=chrlen(c);
  len=max(len,max(R.end(R.chr==c)));
  r = false(len,1);
  for i=1:slength(R),if R.chr(i)==c
    r(R.start(i):R.end(i))=true;
  end, end
  outstem = [outd '/chr' num2str(c)];
  save([outstem '.mat'],'r');
  out = fopen([outstem '.txt'],'wt');
  fprintf(out,'%d\n',r);
  fclose(out);
end


%%%%   (defined by 29-mammal alignment from Gaddy)

load('/xchip/cga1/lawrence/mm/analysis/20100318_cons/regulatory_regions_29mammal.mat','R');

chrlen = load_chrlen;
outd = '/xchip/cga1/lawrence/db/regulatory29';
mkdir(outd);

for c=1:24,disp(c);
  len=chrlen(c);
  idx = find(R.chr==c);
  if ~isempty(idx), len=max(len,max(R.end(idx))); end
  r = false(len,1);
  for i=1:slength(R),if R.chr(i)==c
    r(R.start(i):R.end(i))=true;
  end, end
  outstem = [outd '/chr' num2str(c)];
  save([outstem '.mat'],'r');
  out = fopen([outstem '.txt'],'wt');
  fprintf(out,'%d\n',r);
  fclose(out);
end








