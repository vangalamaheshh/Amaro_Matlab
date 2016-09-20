function choose_GC_bias_matched_lanes_for_segseq(inmatfile,laneinfofile,outmatfile)
% choose_GC_bias_matched_lanes_for_segseq(inmatfile,laneinfofile,outmatfile)
%
% Mike Lawrence 2010-06

fprintf('Loading reads\n');
L = load_struct(laneinfofile);
fc = load(inmatfile);

fprintf('Calculating GC bias profiles\n');
gct = -ones(size(fc.READT,1),1);
gcn = -ones(size(fc.READN,1),1);
for chr=1:23, fprintf('chr%d ',chr);
  load(['/xchip/cga1/lawrence/db/gc100/chr' num2str(chr) '.mat'],'gc');
  idx = find(fc.READT.chr==chr); pos = min(length(gc),fc.READT.pos1(idx)+50); gct(idx) = gc(pos);
  idx = find(fc.READN.chr==chr); pos = min(length(gc),fc.READN.pos1(idx)+50); gcn(idx) = gc(pos);
end, fprintf('\n');

first_tlane = double(min(fc.READT.lane)); last_tlane = double(max(fc.READT.lane));
first_nlane = double(min(fc.READN.lane)); last_nlane = double(max(fc.READN.lane));
n = hist2d_fast(double(gcn),double(fc.READN.lane),0,100,first_nlane,last_nlane);
t = hist2d_fast(double(gct),double(fc.READT.lane),0,100,first_tlane,last_tlane);
tn = [t n];
tlane = [(first_tlane:last_tlane)';nan(last_nlane-first_nlane+1,1)];
nlane = [nan(last_tlane-first_tlane+1,1);(first_nlane:last_nlane)'];

tn2 = bsxfun(@rdivide,tn,sum(tn,1));

fprintf('Clustering lanes according to GC bias\n');
done = false;

maxclust = 10;
while(~done)
  % perform clustering
  fprintf('Clustering with k=%d\n',maxclust);
  c = cluster(linkage(pdist(tn2','correlation')),'maxclust',maxclust);
  % draw figure
  [tmp ord] = sort(c);imagesc(tn2(:,ord)); colorbar; ylim([20 65]);
  idx = find(diff(c(ord))); for i=1:length(idx), line(idx(i)+[0.5 0.5],[0.5 100.5],'color',[0 0 0]); end
  set(gca,'xtick',1:length(c),'xticklabel',nansub({'T','N'},isnorm(ord)+1));

  % choose which cluster
  % (should be large, should be nearly half+half T+N reads)
  lanereads = sum(tn,1);
  istum = ~isnan(tlane(1:length(c)));
  isnorm = ~isnan(nlane(1:length(c)));
  nclust = max(c);
  for i=1:nclust
    ctreads(i,1) = sum(lanereads(istum & c==i));
    cnreads(i,1) = sum(lanereads(isnorm & c==i));
  end
  prod = ctreads.*cnreads;
  if all(prod==0)
    fprintf('No cluster contains both tumor and normal reads\n');
    maxclust = maxclust - 1;
  else
    [tmp ord] = sort(prod,'descend');
    clusterno = ord(1);
    done = true;
    fprintf('Chose cluster %d\n',clusterno);
  end
end

fprintf('Type "return" to continue\n');
keyboard

% save output file
fprintf('Saving output file\n');
idx = find(c==clusterno);
tidx = find(ismember(fc.READT.lane,tlane(idx)));
nidx = find(ismember(fc.READN.lane,nlane(idx)));

READT = [];
READT.chr = int8(fc.READT.chr(tidx));
READT.pos = int32(fc.READT.pos1(tidx));
READT.lane = int16(fc.READT.lane(tidx));

READN = [];
READN.chr = int8(fc.READN.chr(nidx));
READN.pos = int32(fc.READN.pos1(nidx));
READN.lane = int16(fc.READN.lane(nidx));

load('/home/radon01/lawrence/CancerGenomeAnalysis/trunk/segseq/HG18_N36_D2_WINDOWS_100K');
%addpath /home/radon01/lawrence/CancerGenomeAnalysis/trunk/segseq
RATIOS = calc_ratios_from_reads(READN,READT,WINDOWS);

save(outmatfile,'READT','READN','RATIOS','-v7.3');


























return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% other approaches











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% (other early code, not used)
plot(tn2([1:50 52:end],:));

plot(cumsum(tn2([1:50 52:end],:)));

imagesc(tn2);
set(gca,'xtick',1:8,'xticklabel',{'T-sm','T-sm','T-lg','T-lg','N-sm','N-sm','N-lg','N-lg'});

[val idx] = max(tn2);

laneidx=  find(sum(tn,1)>1000);
d = 1-dist(tn2(:,laneidx)',[],'correlation');

f = nan(1,size(tn2,2)); for i=1:size(tn2,2), idx = find(tn2(:,i)>0.03,1); if ~isempty(idx), f(i) = idx; end, end
g = nan(1,size(tn2,2)); for i=1:size(tn2,2), idx = find(tn2(:,i)>0.035,1); if ~isempty(idx), g(i) = idx; end, end
[tmp ord] = sortrows([f' g']);imagesc(tn2(:,ord)); colorbar; ylim([20 55]); png

c = cumsum(tn2);
h = nan(1,size(tn2,2)); for i=1:size(tn2,2), idx = find(c(:,i)>=0.1,1); if ~isempty(idx), h(i) = idx; end, end
j = nan(1,size(tn2,2)); for i=1:size(tn2,2), idx = find(c(:,i)>=0.3,1); if ~isempty(idx), j(i) = idx; end, end
k = nan(1,size(tn2,2)); for i=1:size(tn2,2), idx = find(c(:,i)>=0.5,1); if ~isempty(idx), k(i) = idx; end, end
y = h+j+k;
[tmp ord] = sort(y);imagesc(tn2(:,ord)); colorbar; ylim([20 55]); png


