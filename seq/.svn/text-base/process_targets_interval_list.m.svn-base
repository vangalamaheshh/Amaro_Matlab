function process_targets_interval_list(infile,outfile)
% process_targets_interval_list(infile,outfile)
%
% Mike Lawrence 2009

build = 'hg18'

T = load_target_interval_list(infile)

% find out what genes they are

R = load_refseq(build);
R.chr = convert_chr(R.chr);

nt = slength(T);
T.gene = cell(nt,1);

% pass 1: if hits only one gene's exon(s)
for i=1:nt
  idx = find(R.chr==T.chr(i) & R.code_start <= T.end(i) & R.code_end >= T.start(i));
  g = unique(R.gene(idx));
  if isempty(g)  fprintf('%d. None found\n',i);
  elseif length(g)>1
    exonic = []; for j=1:length(idx), for e=1:R.n_exons(idx(j))
        if R.exon_starts{idx(j)}(e) <= T.end(i) & R.exon_ends{idx(j)}(e) >= T.start(i)
           exonic = [exonic; idx(j)];
    end,end,end
    g_exonic = unique(R.gene(exonic));
    if isempty(g_exonic)  fprintf('%d. None exonic found\n',i);
    elseif length(g_exonic)>1, fprintf('%d. Multiple exonic found:',i);
      for j=1:length(g_exonic), fprintf(' %s',g_exonic{j}); end, fprintf('\n');
    else T.gene{i}=g_exonic{1}; end
  else T.gene{i}=g{1}; end
end
fprintf('%d targets with unknown gene\n',sum(cellfun('isempty',T.gene)));

% pass 2: annotate pseudo-mutations
len = T.end-T.start+1;
for pass=1:3
  idx = find(cellfun('isempty',T.gene));
  tmp = [];
  tmp.chr = T.chr(idx);
  switch(pass)
    case 1, tmp.pos = round((T.start(idx)+T.end(idx))/2);
    case 2, tmp.pos = T.start(idx)+round(0.1*len(idx));
    case 3, tmp.pos = T.start(idx)+round(0.9*len(idx));
  end
  tmp.newbase = cell(length(idx),1);
  mut=[];
  mut.from = {'A';'C';'G';'T'};
  mut.to = {'C';'G';'T';'A'};
  for c=1:24
    jdx = find(tmp.chr==c);
    tmp.newbase(jdx) = mapacross(upper(genome_region(c,tmp.pos(jdx))),mut.from,mut.to);
  end
  tmp = classify_muts(tmp,struct('impute_promoters',false),R);
  jdx = grep('miRNA|Missense|Synonymous|Nonsense|Splice|Read-through',tmp.type,1);
  T.gene(idx(jdx)) = tmp.gene(jdx); 
  fprintf('%d targets with unknown gene\n',sum(cellfun('isempty',T.gene)));
end

% pass 3: if hits only one gene's transcribed area
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if isempty(g)  fprintf('%d. None found\n',i);
  elseif length(g)>1, fprintf('%d. Multiple found\n',i); disp(g);
  else fprintf('%d. Matched: %s\n',i,g{1}); T.gene{i}=g{1}; end
end
fprintf('%d targets with unknown gene\n',sum(cellfun('isempty',T.gene)));

% pass 4: SNOR's get priority (because they're so small)
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  sn = grep('^SNOR',g);
  if ~isempty(sn), fprintf('%d. Matched: %s\n',i,sn{1}); T.gene{i}=sn{1}; end
end
fprintf('%d targets with unknown gene\n',sum(cellfun('isempty',T.gene)));

% pass 5: first gene on list
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if ~isempty(g), fprintf('%d. picking %s\n',i,g{1}); T.gene{i}=g{1}; end
end
fprintf('%d targets with unknown gene\n',sum(cellfun('isempty',T.gene)));

% rest are "Unknown"
T.gene = fillblanks(T.gene,'Unknown');

% add GC content
nt = slength(T);
T.gc = nan(nt,1);
for chr=1:24
  fprintf('Chr %d\n',chr);
  ref = upper(genome_region(chr,1,inf));
  len = length(ref);
  gc = nan(len,1);
  gc(ref=='A' | ref=='T') = 0;
  gc(ref=='C' | ref=='G') = 1;
  idx = find(T.chr==chr);
  for j=1:length(idx), i=idx(j);
    st = max(1,T.start(i));
    en = min(len,T.end(i));
    T.gc(i) = nanmean(gc(st:en));
  end
end

% save
T = keep_fields(T,{'gene','chr','start','end','gc'});
save_struct(T,outfile,'no_headers');



return

%%%%%% OLD METHODS



% pass 2: if hits exons of adjacent gene on list
while(1)
n = sum(cellfun('isempty',T.gene));
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if length(g)>1
    adjacent = {};
    if i>1, adjacent = [adjacent; T.gene{i-1}]; end
    if i<nt, adjacent = [adjacent; T.gene{i+1}]; end
    g_adjacent = intersect(g,adjacent);
    if isempty(g_adjacent)  fprintf('%d. None adjacent found\n',i);
    elseif length(g_adjacent)>1, fprintf('%d. Multiple adjacent found:',i);
      for j=1:length(g_adjacent), fprintf(' %s',g_adjacent{j}); end, fprintf('\n');
    else T.gene{i}=g_adjacent{1}; end
  end
end
n_new = sum(cellfun('isempty',T.gene));
if (n_new==n) break; end
end

% pass 3: pick first overlapping gene
for i=1:nt
  if ~isempty(T.gene{i}), continue; end
  idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i) & R.tx_end >= T.start(i));
  g = unique(R.gene(idx));
  if length(g)>1, T.gene{i}=g{1}; end
end


% pass 4: nearby genes
for margin=[1000,2000,5000,10000,20000,50000,100000,200000,500000,10000000]
  for i=1:nt
    if ~isempty(T.gene{i}), continue; end
    idx = find(R.chr==T.chr(i) & R.tx_start <= T.end(i)+margin & R.tx_end >= T.start(i)-margin);
    g = unique(R.gene(idx));
    if length(g)>1, T.gene{i}=g{1}; end
  end
end

sum(cellfun('isempty',T.gene))  % 0
T.gene = fillblanks(T.gene,'Unknown');

% add GC content

nt = slength(T);
T.gc = nan(nt,1);
for chr=1:24
  fprintf('Chr %d\n',chr);
  ref = upper(genome_region(chr,1,inf));
  len = length(ref);
  gc = nan(len,1);
  gc(ref=='A' | ref=='T') = 0;
  gc(ref=='C' | ref=='G') = 1;
  idx = find(T.chr==chr);
  for j=1:length(idx), i=idx(j);
    st = max(1,T.start(i));
    en = min(len,T.end(i));
    T.gc(i) = nanmean(gc(st:en));
  end
end

% save
T = keep_fields(T,{'gene','chr','start','end','gc'});
save_struct(T,outfile,'no_headers');
