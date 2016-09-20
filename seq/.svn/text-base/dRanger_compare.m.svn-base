function Y = dRanger_compare(sample_list,tolerance)
% OBSOLETE: replaced by dRanger_intersect

if ~exist('tolerance','var'), tolerance = 100; end

ns = length(sample_list);
X = cell(ns,1);
B = [];
for i=1:ns
  sample = sample_list{i};
  fprintf('Loading sample %s\n',sample);
  direc = ['/xchip/tcga_scratch/lawrence/' sample];
  name2 = upper(regexprep(sample,'/','-'));
  fname = [direc '/' name2 '_dRanger_results_bfiltered.txt'];
  if ~exist(fname,'file'), error('Can''t find file %s',fname);end
  X{i} = load_struct(fname);
  X{i} = make_numeric(X{i},{'filterHCL','filterB'});
  X{i} = reorder_struct(X{i},~X{i}.filterHCL & ~X{i}.filterB);
  X{i} = make_numeric(X{i},{'num','chr1','chr2','pos1','pos2','str1','str2','tumreads','normreads'});
  B = [B; [repmat(i,slength(X{i}),1) X{i}.num X{i}.chr1 X{i}.pos1 X{i}.str1 X{i}.chr2 X{i}.pos2 X{i}.str2]];
end

B1 = B;
B1(:,[4 7]) = round(B1(:,[4 7]) / tolerance);
length(B1)
length(unique(B1(:,3:end),'rows'))

[u ui uj] = unique(B1(:,3:8),'rows');
ny = length(u);
Y = [];
Y.chr1=B1(ui,3); Y.str1=B1(ui,5);
Y.chr2=B1(ui,6); Y.str2=B1(ui,8);

Y.idx = nan(ny,ns); Y.num = nan(ny,ns);
Y.ntum = nan(ny,ns); Y.nnorm = nan(ny,ns);
Y.pos1 = nan(ny,ns); Y.pos2 = nan(ny,ns);
for i=1:ny
  j = find(uj==i);
  for s=1:ns
    k = j(find(B1(j,1)==s,1));
    if ~isempty(k)
      Y.num(i,s) = B(k,2); Y.idx(i,s) = find(X{s}.num==Y.num(i,s));
      Y.ntum(i,s) = X{s}.tumreads(Y.idx(i,s)); Y.nnorm(i,s) = X{s}.normreads(Y.idx(i,s));
      Y.pos1(i,s)=X{s}.pos1(Y.idx(i,s)); Y.pos2(i,s)=X{s}.pos2(Y.idx(i,s));
    end
  end
end

Y.dpos1 = Y.pos1-repmat(round(nanmean(Y.pos1,2)),1,ns);
Y.dpos2 = Y.pos2-repmat(round(nanmean(Y.pos2,2)),1,ns);

Y.nsamps = sum(~isnan(Y.idx),2);

sortrows([Y.nsamps Y.ntum Y.nnorm Y.dpos1 Y.dpos2])

fprintf('%d rearrangements\n%d bp tolerance\n%d unique rearrangements',...
  length(B),tolerance,ny);
count(Y.nsamps)

