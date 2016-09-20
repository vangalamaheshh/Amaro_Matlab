function process_TN_wigs(bd,cen)
% Mike Lawrence 2009-12-11

% bd = '/xchip/tcga/gbm/analysis/lawrence/tcga/gbm/20090629/vcf';
% cen = {'broad','baylor','washu'};

try
  g = Genome;
catch me
  fprintf('Please add file "java.opts" in the startup directory, containing the line "-Xmx5g"\n');
  throw(me)
end

for c=1:length(cen)
  cd = [bd '/' cen{c}]; d = dir([cd '/TCGA-*.wig']); d = list2cell(d.name);
  d = setdiff(d,grep('max|somcov',d));
  q = parse(d,'^(TCGA-..-....)-(.)(.*)\.wig$',{'patient','tn','suffix'});
  [u ui uj] = unique(q.patient);
  for i=1:length(u)
    fprintf('Processing %s sample %s\n',cen{c},u{i});
    tidx = find(uj==i & strcmp(q.tn,'0')); nidx = find(uj==i & strcmp(q.tn,'1'));
    if length(tidx)<1 || length(nidx)<1
      fprintf('  Skipping: %d tumor and %d normal files\n',length(tidx),length(nidx)); continue;
    end
    % if more than one tumor file, unite them (or)
    if length(tidx)>1
      fprintf('  OR-ing tumor files:\n');
      g.clear();
      for j=1:length(tidx)
        fname = d{tidx(j)};
        g.loadWiggle([cd '/' fname],'or');
        fprintf('    File (%d) %s   \ttotContents=%ld\n',j,fname,g.totContents);
      end
      tname = [u{i} '-tumor.wig'];
      g.saveWiggle([cd '/' tname]);
    else
      tname = d{tidx};
    end
    % if more than one normal file, unite them (or)
    if length(nidx)>1
      fprintf('  OR-ing normal files:\n');
      g.clear();
      for j=1:length(tidx)
        fname = d{nidx(j)};
        g.loadWiggle([cd '/' fname],'or');
        fprintf('    File (%d) %s   \ttotContents=%ld\n',j,fname,g.totContents);
      end
      nname = [u{i} '-normal.wig'];
      g.saveWiggle([cd '/' nname]);
    else
      nname = d{nidx};
    end
    % unite tumor+normal (and)
    fprintf('  AND-ing tumor and normal file:\n');
    g.clear();
    g.loadWiggle([cd '/' tname]);
    fprintf('    Tumor file %s   \ttotContents=%ld\n',tname,g.totContents);
    g.loadWiggle([cd '/' nname],'and');
    fprintf('    Normal file %s   \ttotContents=%ld\n',nname,g.totContents);
    pname = [u{i} '-somcov.wig'];
    fprintf('    AND-ed patient file %s   \ttotContents=%ld\n',pname,g.totContents);
    g.saveWiggle([cd '/' pname]);
  end
end

% for each patient, unite three centers' wig files

q = cell(length(cen),1);
for c=1:length(cen)
  cd = [bd '/' cen{c}]; d = dir([cd '/TCGA*somcov.wig']); d = list2cell(d.name);
  q{c} = parse(d,'^(TCGA-..-....)-somcov.wig$',{'patient'});
  q{c}.cen = c*ones(length(d),1);
end
q=concat_structs(q);
[u ui uj] = unique(q.patient);
ud = [bd '/union'];
if ~exist(ud,'dir'), mkdir(ud); end
for i=1:length(u)
  fprintf('Processing sample %s\n',u{i});
  idx = find(uj==i);
  fprintf('  OR-ing files:\n');
  g.clear();
  for j=1:length(idx)
    fname = [cen{q.cen(idx(j))} '/' u{i} '-somcov.wig'];
    g.loadWiggle([bd '/' fname],'or');
    fprintf('    File (%d) %s   \ttotContents=%ld\n',j,fname,g.totContents);
  end
  uname = [u{i} '-somcov.wig'];
  fprintf('    OR-ed file %s   \ttotContents=%ld\n',uname,g.totContents);
  g.saveWiggle([ud '/' uname]);
end

