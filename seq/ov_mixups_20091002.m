%%
%%  DETECTING MIXUPS IN CAPTURE DATA
%%


%%%%  2009-09-28
% mixup detection using hybrid selection metrics

% GET SAMPLE LIST

%%% RUN THE FOLLOWING COMMAND IN ISOLATION TO AVOID STDIN->STDOUT INTERFERENCE  %%%
[tmp list] = system('ls -ld /xchip/cga1/lawrence/ov/*/capture/*.bam.lanetable');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

list = setdiff(split(list,char(10)),{''});
fname = regexprep(list,'^.*(/xchip/.*lanetable$)','$1');

% GET LANETABLES

L = cell(length(fname),1);
for i=1:length(fname), L{i} = load_struct(fname{i}); end
L = concat_structs(L);
tmp = parse(L.SM,'^TCGA-\d\d-(\d\d\d\d)-(\d)\d.-...$',{'indiv','is_norm'});
L = merge_structs({L,tmp});
tmp = parse(L.PU,'^(.........).*\.(\d)$',{'fc_name_long','fc_lane'});
tmp.fc_name = regexprep(tmp.fc_name_long,'^(.....).*$','$1');
L = merge_structs({L,tmp});

% LOOK UP HYBRID-SELECTION METRICS IN PICARD

metrics = {'BAIT_SET','PCT_SELECTED_BASES','ON_BAIT_VS_SELECTED','FOLD_ENRICHMENT','FOLD_80_BASE_PENALTY','HS_LIBRARY_SIZE'};
fprintf('Loading hybrid-selection metrics: Flowcell ');
[f fi fj] = unique(L.fc_name_long);
M = repmat({''},slength(L),length(metrics));
for i=1:length(f), fprintf('%d/%d ',i,length(f));
  fclanes = find(fj==i);
  dname = ['/seq/picard/' f{i}];
  d = dir([dname '/*_*']);
  d = d(find(cat(1,d.isdir)));
  [tmp ord] = sort(cat(1,d.datenum),'descend');
  for j=1:length(fclanes), k=fclanes(j);
    found = false;
    for di=1:length(d)   % try the subdirectories in order of newest->oldest
      fname = [dname '/' d(ord(di)).name '/' L.fc_lane{k} '/' L.LB{k} '/' f{i} '.' L.fc_lane{k} '.hybrid_selection_metrics'];
      if exist(fname,'file')
        H = load_picard_struct(fname);
        if slength(H)>1, error('too many lines in %s',fname); end
        for m=1:length(metrics)
          if isfield(H,metrics{m}), M(k,m) = getfield(H,metrics{m}); end
        end
        found = true;
        break;
      end
    end
    if ~found, fprintf('\nNo hybrid_selection_metrics found for flowcell %s lane %s\n',f{i},L.fc_lane{k}); end
  end
end, fprintf('\n');
for m=1:length(metrics), L = setfield(L,metrics{m},M(:,m)); end

tn = {'t','n'};
for i=1:slength(L), L.sample{i,1} = [L.indiv{i} tn{str2double(L.is_norm{i})+1}]; end

save('/xchip/tcga_scratch/lawrence/ov/analysis/20090929/L.mat','L');
load('/xchip/tcga_scratch/lawrence/ov/analysis/20090929/L.mat');

% EXPLORE METRICS

[LB.name lbi lbj] = unique(L.LB);
LB.SM = L.SM(lbi);
LB.indiv = L.indiv(lbi);
LB.is_norm = L.is_norm(lbi);
LB.sample = L.sample(lbi);
for i=1:slength(LB)
  idx = find(lbj==i);
  LB.nlanes(i,1) = length(idx);
  LB.psb_std(i,1) = nanstd(str2double(L.PCT_SELECTED_BASES(idx)));
  LB.obs_std(i,1) = nanstd(str2double(L.ON_BAIT_VS_SELECTED(idx)));
  LB.fe_std(i,1) = nanstd(str2double(L.FOLD_ENRICHMENT(idx)));
  LB.f80bp_std(i,1) = nanstd(str2double(L.FOLD_80_BASE_PENALTY(idx)));
  LB.hls_std(i,1) = nanstd(str2double(L.HS_LIBRARY_SIZE(idx)));
end

%% BLACKLIST

BL=[];
BL.fcl = load_lines('/xchip/cga1/gdac-prod/genepattern/jobResults/17447/lane_black_list.txt');
for i=1:slength(L), L.fcl{i,1} = [L.fc_name{i} '.' L.fc_lane{i}]; end
BL.sample = map_across(BL.fcl,L.fcl,L.sample);
BL.library = map_across(BL.fcl,L.fcl,L.LB);
bl = unique(BL.sample);
fprintf('%s\n',bl{:});

col = zeros(slength(LB),3);   %  black for tumors
col (find(strcmp(LB.is_norm,'1')),:) = 0.6;   % grey for normals
col(ismember(LB.name,BL.library),1) = 1;  % red

% obs
t = metrics{3};
x = LB.nlanes;
y = LB.obs_std;
x = x + (0.3*rand(slength(LB),1)-0.15);
scatter(x,y,[],col)
bad = find(y>0.1);
for j=1:length(bad),i=bad(j);
  text(x(i)-0.55,y(i),LB.sample{i},'color',col(i,:));
end
xlabel('number of lanes from library','fontsize',15);
ylabel('stdev','fontsize',15);
title(regexprep(t,'_','\\_'),'fontsize',15);
set(gca,'xtick',1:6);
line([1.3 6.5],[.12 .05],'color',[1 0 0]);

% obs: actual data
good4 = setdiff(find(LB.nlanes==4),bad);
good6 = setdiff(find(LB.nlanes==6),bad);
good = [good4(1:4);good6(1:4)];
show = [good;bad];
x = []; y = []; c = []; name = {}; idx=1;
for j=1:length(show), i=show(j);
  lidx = find(strcmp(L.LB,LB.name{i}));
  name{idx} = LB.sample{i};
  c = [c;repmat(col(i,:),length(lidx),1)];
  x = [x;repmat(idx,length(lidx),1)];
  y = [y;L.ON_BAIT_VS_SELECTED(lidx)];
  idx=idx+1;
end
if ~isnumeric(y), y = str2double(y); end
scatter(x,y,[],c);
xlabel('library','fontsize',15);
ylabel('ON\_BAIT\_VS\_SELECTED','fontsize',15);
set(gca,'xtick',1:max(x)+1);
ylim([0 1.1]);xlim([0 max(x)+0.5]);
for i=1:max(x)
  line([i i],[min(y(x==i)) max(y(x==i))]);
  text(i,1,name{i},'rotation',90,'color',c(find(x==i,1),:));
end



% psb
t = metrics{2};
x = LB.nlanes;
y = LB.psb_std;
x = x + (0.3*rand(slength(LB),1)-0.15);
scatter(x,y,[],col)
bad = find(y>0.14 | (x>2.5 & y>0.1));
for j=1:length(bad),i=bad(j);
  text(x(i)-0.55,y(i),LB.sample{i},'color',col(i,:));
end
xlabel('number of lanes from library','fontsize',15);
ylabel('stdev','fontsize',15);
title(regexprep(t,'_','\\_'),'fontsize',15);
set(gca,'xtick',1:6);
line([1.3 6.5],[.14 .06],'color',[1 0 0]);

% psb: actual data
good4 = setdiff(find(LB.nlanes==4),bad);
good6 = setdiff(find(LB.nlanes==6),bad);
good = [good4(1:4);good6(1:4)];
show = [good;bad];
x = []; y = []; c = []; name = {}; idx=1;
for j=1:length(show), i=show(j);
  lidx = find(strcmp(L.LB,LB.name{i}));
  name{idx} = LB.sample{i};
  c = [c;repmat(col(i,:),length(lidx),1)];
  x = [x;repmat(idx,length(lidx),1)];
  y = [y;L.PCT_SELECTED_BASES(lidx)];
  idx=idx+1;
end
if ~isnumeric(y), y = str2double(y); end
scatter(x,y,[],c);
xlabel('library','fontsize',15);
ylabel('PCT\_SELECTED\_BASES','fontsize',15);
set(gca,'xtick',1:max(x)+1);
ylim([0 1.1]);xlim([0 max(x)+0.5]);
for i=1:max(x)
  line([i i],[min(y(x==i)) max(y(x==i))]);
  text(i,1,name{i},'rotation',90,'color',c(find(x==i,1),:));
end



% fe
t = metrics{4};
x = LB.nlanes;
y = LB.fe_std;
x = x + (0.3*rand(slength(LB),1)-0.15);
scatter(x,y,[],col)
bad = find(y>90 | (x>3.5 & y>20));
for j=1:length(bad),i=bad(j);
  xoff=0; yoff=0;
  if strcmp(LB.sample{i},'1328n'), yoff=+1; end
  if ismember(LB.sample{i},{'0920t','0979t'}), yoff=-1; end
  if strcmp(LB.sample{i},'1321n'), xoff=+0.6; end
  text(x(i)-0.55+xoff,y(i)+yoff,LB.sample{i},'color',col(i,:));
end
xlabel('number of lanes from library','fontsize',15);
ylabel('stdev','fontsize',15);
title(regexprep(t,'_','\\_'),'fontsize',15);
set(gca,'xtick',1:6);
line([1.3 6.5],[77 7],'color',[1 0 0]);




% fe: actual data
good4 = setdiff(find(LB.nlanes==4),bad);
good6 = setdiff(find(LB.nlanes==6),bad);
good = [good4(1:4);good6(1:4)];
show = [good;bad];
x = []; y = []; c = []; name = {}; idx=1;
for j=1:length(show), i=show(j);
  lidx = find(strcmp(L.LB,LB.name{i}));
  name{idx} = LB.sample{i};
  c = [c;repmat(col(i,:),length(lidx),1)];
  x = [x;repmat(idx,length(lidx),1)];
  y = [y;L.FOLD_ENRICHMENT(lidx)];
  idx=idx+1;
end
if ~isnumeric(y), y = str2double(y); end
scatter(x,y,[],c);
xlabel('library','fontsize',15);
ylabel('FOLD\_ENRICHMENT','fontsize',15);
set(gca,'xtick',1:max(x)+1);
ylim([0 350]);xlim([0 max(x)+0.5]);
for i=1:max(x)
  line([i i],[min(y(x==i)) max(y(x==i))]);
  ypos = 140; if i==9, ypos = 130; end
  text(i,ypos,name{i},'rotation',90,'color',c(find(x==i,1),:));
end










% f80bp
t = metrics{5};
x = LB.nlanes;
y = LB.f80bp_std;
x = x + (0.3*rand(slength(LB),1)-0.15);
scatter(x,y,[],col)
xlabel('number of lanes from library','fontsize',15);
ylabel('stdev','fontsize',15);
title(regexprep(t,'_','\\_'),'fontsize',15);
set(gca,'xtick',1:6);


% hls
t = metrics{6};
x = LB.nlanes;
y = log10(LB.hls_std);
x = x + (0.3*rand(slength(LB),1)-0.15);
scatter(x,y,[],col)
xlabel('number of lanes from library','fontsize',15);
ylabel('log(stdev)','fontsize',15);
title(regexprep(t,'_','\\_'),'fontsize',15);
set(gca,'xtick',1:6);





%%%%%%  2009-10-01
%%%%%%
%%%%%%  CN QC using SNP-array ground truth

%%  get sample list

samples = build_OV_capture_list;

% prepare absolute coverage in each whole-exome target
exome_specs = ['/xchip/tcga_scratch/lawrence/capture/'...
  'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP_noheader.txt'];
cov_filestem = 'exome_cov_per_lane';
RegionCovPerLane(samples,exome_specs,cov_filestem);

% load lanetables and regioncov
R=[]; L={}; Cr={}; idx=0;
tn={'tumor','normal'};
for i=1:length(samples),fprintf('%d/%d ',i,length(samples));
  for j=1:2, idx=idx+1;
    L{idx} = load_struct(['/xchip/tcga_scratch/lawrence/' samples{i} '/' tn{j} '.bam.lanetable']);
    L{idx}.sample = repmat(samples(i),slength(L{idx}),1);
    L{idx}.istum = repmat(strcmp(tn{j},'tumor'),slength(L{idx}),1);
    x = read_table(['/xchip/tcga_scratch/lawrence/' samples{i} '/exome_cov_per_lane_' tn{j} '.txt'],...
      ['%s' repmat('%f',1,3+slength(L{idx}))],char(9),0);
    if isempty(R), R.gene = x.dat{1};R.chr = x.dat{2};R.start = x.dat{3};R.end = x.dat{4}; end
    Cr{idx} = cat(2,x.dat{5:end});
  end
end,fprintf('\n');
L = concat_structs(L);
Cr = cat(2,Cr{:});
nl = slength(L);
% sort regions by position
[R ord] = sort_struct(R,{'chr','start'});
Cr = Cr(ord,:);

% get other region info
R = load_struct(['/xchip/tcga_scratch/lawrence/capture/'...
  'whole_exome_refseq_coding.targets.interval_list.GENESGC.WITHMEMBERSHIP.txt'],...
  '%s%f%f%f%f%f%f');
R = sort_struct(R,{'chr','start'});

gc_min = 0.45;
gc_max = 0.55;
R2 = reorder_struct(R,R.gc>=gc_min & R.gc<=gc_max);

% coverage: collapse regions to genes
G = [];
[G.name gi gj] = unique(R2.gene);
G.chr = R2.chr(gi);
G.membership = R2.membership(gi);
ng = slength(G);
Cg = zeros(ng,nl);
for i=1:ng
  idx = find(gj==i);
  Cg(i,:) = sum(Cr(idx,:),1);
  G.start(i,1) = min(R2.start(idx));
  G.end(i,1) = max(R2.end(idx));
  G.len(i,1) = sum(R2.len(idx));
  G.gc(i,1) = weighted_mean(R2.gc(idx),R2.len(idx));
end
G.mid = (G.start+G.end)/2;
% sort genes by position
[G ord] = sort_struct(G,{'chr','start'});
Cg = Cg(ord,:);

% retrieve SNP array ground truth
% (segmented data)

seg = {};
for batch=9:20
  sd = seg_data();
  try
    fprintf('Trying to retrieve batch %d\n',batch);
    sd = seg_data();
    seg{batch} = sd.get_batch_data('broad',batch);
  catch
    fprintf('  Batch not found.\n');
  end
end
seg = cat(1,seg{~cellfun('isempty',seg)});
Z = [];
Z.sample = seg.Sample;
Z.chr = seg.Chromosome;
Z.start = seg.Start;
Z.end = seg.End;
Z.nprobes = seg.Num_Probes;
Z.segmean = seg.Segment_Mean;

%     Sample          Chromosome    Start         End           Num_Probes    Segment_Mean
% (sorted by Chromsome,Start)


%%%% FILTER AND NORMALIZE DATA

% look only at 6k data
lidx = grep('6k',L.baitset,1);
gidx = find(G.membership>1);
L6 = reorder_struct(L,lidx);
G6 = reorder_struct(G,gidx);
Cg6 = Cg(gidx,lidx);
ltot = sum(Cg6,1);
gtot = sum(Cg6,2);
% look only at lanes with >1e6 reads
lidx = find(ltot>1e6);
L6 = reorder_struct(L6,lidx);
Cg6 = Cg6(:,lidx);
ltot = sum(Cg6,1);
% normalize all lanes to 1e6 reads
Cg6 = round(bsxfun(@rdivide,1e6*Cg6,ltot));
% normalize to gene length
Cg6 = round(bsxfun(@rdivide,100*Cg6,G6.len));
ltot = sum(Cg6,1);
gtot = sum(Cg6,2);
% look at only genes with typical coverage-per-length
gidx = find(gtot>1000 & gtot<10000);
G6 = reorder_struct(G6,gidx);
Cg6 = Cg6(gidx,:);
ltot = sum(Cg6,1);
gtot = sum(Cg6,2);
% get rid of weird lane #149
lidx = setdiff(1:slength(L6),149);
L6 = reorder_struct(L6,lidx);
Cg6 = Cg6(:,lidx);
% sort lanes
[L6 ord] = sort_struct(L6,{'istum','sample'},[-1 1]);
Cg6 = Cg6(:,ord);
% divide by average of normals
a = mean(Cg6(:,~L6.istum),2);
Cg6 = bsxfun(@rdivide,Cg6,a);


save('/xchip/tcga_scratch/lawrence/ov/analysis/20090929/20091002_cn_qc_ALL.mat');
load('/xchip/tcga_scratch/lawrence/ov/analysis/20090929/20091002_cn_qc_ALL.mat');


%%%  reformat segments to genes

seg.shortname = regexprep(seg.Sample,'^(TCGA-..-....-..).*$','$1');
S = [];
[S.name si sj] = unique(seg.shortname);
ns = slength(S);
ng6 = slength(G6);
A6 = nan(ng6,ns);
for s=1:ns, if ~mod(s,10), fprintf('%d/%d ',s,ns); end
  sdx = find(sj==s);
  gidx = find(G6.chr==Z.chr(1) & G6.start>=Z.start(1),1);
  for j=1:length(sdx), i=sdx(j);
    if Z.nprobes(i)<500, continue; end
    n = find(G6.chr(gidx:end)>Z.chr(i) | G6.mid(gidx:end)>Z.end(i),1);
    if ~isempty(n), A6(gidx:gidx+n-1,s) = Z.segmean(i); gidx=gidx+n-1;
    else A6(gidx:end,s) = Z.segmean(i);
end,end,end,fprintf('\n');



% CHROMOSOME CN FIGURES

cen = load_cen;
figure(1);

sample = '1500';
chr=8;
%for chr=1:24
gidx = find(G6.chr==chr & G6.membership>1);
x1 = G6.start(gidx);x2 = G6.end(gidx);
clf,hold on
% genes: tumor
tlidx = find(strcmp(['ov/' sample '/capture'],L6.sample) & L6.istum);
tlaneno = 1;
tlaneno = min(tlaneno,length(tlidx));
y = log(Cg6(gidx,tlidx(tlaneno)));scatter(x1/1e6,y,[],'r','filled');
% genes: normal
nlaneno = 2;
nlaneno = min(nlaneno,length(tlidx));
nlidx = find(strcmp(['ov/' sample '/capture'],L6.sample) & ~L6.istum);
y = log(Cg6(gidx,nlidx(nlaneno)));scatter(x1/1e6,y,[],'k','filled');
% segments from tumor SNP array
sidx = grep(['TCGA-..-' sample '-01'],seg.Sample,1);
q=seg(sidx(seg.Chromosome(sidx)==chr & seg.Num_Probes(sidx)>=500,:),:);
title([seg.Sample{sidx(1)}(1:12) ' tlaneno ' num2str(tlaneno) ' nlaneno ' num2str(nlaneno) ' chr' num2str(chr)]);
line([q.Start q.End]'/1e6,[q.Segment_Mean q.Segment_Mean]','color',[0 0 1],'linewidth',4)
line(repmat(cen(chr,:),2,1)/1e6,[-1.5 -1.5;1.5 1.5],'color','b');line([1 max(x2)]/1e6,[0 0],'color','k');
% segments mapped to genes
aidx = grep(['TCGA-..-' sample '-01'],S.name,1);
scatter(x1/1e6,A6(gidx,aidx),4,'w');
%% graph finished
ylim([-3 3]),hold off
keyboard,end


%%%% RANK-CORRELATIONS



