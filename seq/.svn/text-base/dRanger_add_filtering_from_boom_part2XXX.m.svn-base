function X = dRanger_add_filtering_from_boom_part2(X,t_boomdir,n_boomdir,P)
% Mike Lawrence 2009-09

if ~exist('P','var'), P=[]; end

P=impose_default_value(P,'skip_end2',false);
P=impose_default_value(P,'add_coverage',true);
P=impose_default_value(P,'add_fmapqz',true);
P=impose_default_value(P,'add_avgnmm',true);
P=impose_default_value(P,'add_nuwp',true);

nx = slength(X);
z = nan(nx,1);

tt = {};

if P.add_coverage
  X.coverageT1 = z;
  X.coverageN1 = z;
  if ~P.skip_end2
    X.coverageT2 = z;
    X.coverageN2 = z;
  end
  tt = [tt;'coverage'];
end

if P.add_fmapqz
  X.fmapqzT1 = z;
  X.fmapqzN1 = z;
  if ~P.skip_end2
    X.fmapqzT2 = z;
    X.fmapqzN2 = z;
  end
  tt = [tt;'fmapqz'];
end

if P.add_avgnmm
  X.avgnmmT1 = z;
  X.avgnmmN1 = z;
  if ~P.skip_end2
    X.avgnmmT2 = z;
    X.avgnmmN2 = z;
  end
  tt = [tt;'avgnmm'];
end

if P.add_nuwp
  X.nuwpT1 = z;
  X.nuwpN1 = z;
  if ~P.skip_end2
    X.nuwpT2 = z;
    X.nuwpN2 = z;
  end
  tt = [tt;'nuwp'];
end

if isempty(tt) 
  fprintf('No filtering metrics requested\n');
  return;
end

% first check to make sure no files are missing
%for c=1:24
dir1=dir([t_boomdir '/chr*.nuwp']);
NC=length(dir1);
for c=1:NC
  tstem = [t_boomdir '/chr' num2str(c)];
  nstem = [n_boomdir '/chr' num2str(c)];
  for tti=1:length(tt)
    d = dir([tstem '.' tt{tti}]);
    if isempty(d), error('chr%d.%s missing from tumor dir',c,tt{tti}); end
    d = dir([nstem '.' tt{tti}]);
    if isempty(d), error('chr%d.%s missing from normal dir',c,tt{tti}); end
  end
end

% load and assign metrics

X = make_numeric(X,{'chr1','min1','max1'});
if ~P.skip_end2
  X = make_numeric(X,{'chr2','min2','max2'});
end

for c=1:NC
  fprintf('chr%d\n',c);

  tstem = [t_boomdir '/chr' num2str(c)];
  nstem = [n_boomdir '/chr' num2str(c)];

  d = dir([tstem '.' tt{1}]);
  tlen = d.bytes;
  d = dir([nstem '.' tt{1}]);
  nlen = d.bytes;

  idx1 = find(X.chr1==c);
  st = X.min1(idx1);
  en = X.max1(idx1);
  tst = min(tlen,st);
  ten = min(tlen,en);
  nst = min(nlen,st);
  nen = min(nlen,en);
  fprintf('  end1: loading\n');
  if P.add_coverage
    tc = get_block([tstem '.coverage'],'byte',tst-1,ten-1);
    nc = get_block([nstem '.coverage'],'byte',nst-1,nen-1);
  end
  if P.add_fmapqz
    tz = get_block([tstem '.fmapqz'],'byte',tst-1,ten-1)/100;
    nz = get_block([nstem '.fmapqz'],'byte',nst-1,nen-1)/100;
  end
  if P.add_avgnmm
    tm = get_block([tstem '.avgnmm'],'byte',tst-1,ten-1)/10;
    nm = get_block([nstem '.avgnmm'],'byte',nst-1,nen-1)/10;
  end
  if P.add_nuwp
    tw = get_block([tstem '.nuwp'],'byte',tst-1,ten-1);
    nw = get_block([nstem '.nuwp'],'byte',nst-1,nen-1);
  end
  fprintf('  end1: assigning\n');
  tpos1 = 1;
  npos1 = 1;
  for j=1:length(idx1), i=idx1(j);
    tpos2 = tpos1+ten(j)-tst(j);
    npos2 = npos1+nen(j)-nst(j);
    if st<=tlen
      if P.add_coverage, X.coverageT1(i) = max(tc(tpos1:tpos2)); end
      if P.add_fmapqz, X.fmapqzT1(i) = max(tz(tpos1:tpos2)); end
      if P.add_avgnmm, X.avgnmmT1(i) = max(tm(tpos1:tpos2)); end
      if P.add_nuwp, X.nuwpT1(i) = max(tw(tpos1:tpos2)); end
    end
    if st<=nlen
      if P.add_coverage, X.coverageN1(i) = max(nc(npos1:npos2)); end
      if P.add_fmapqz, X.fmapqzN1(i) = max(nz(npos1:npos2)); end
      if P.add_avgnmm, X.avgnmmN1(i) = max(nm(npos1:npos2)); end
      if P.add_nuwp, X.nuwpN1(i) = max(nw(npos1:npos2)); end
    end
    tpos1=tpos2+1;
    npos1=npos2+1;
  end

 if ~P.skip_end2
  idx2 = find(X.chr2==c);
  st = X.min2(idx2);
  en = X.max2(idx2);
  tst = min(tlen,st);
  ten = min(tlen,en);
  nst = min(nlen,st);
  nen = min(nlen,en);
  fprintf('  end2: loading\n');
  if P.add_coverage
    tc = get_block([tstem '.coverage'],'byte',tst-1,ten-1);
    nc = get_block([nstem '.coverage'],'byte',nst-1,nen-1);
  end
  if P.add_fmapqz
    tz = get_block([tstem '.fmapqz'],'byte',tst-1,ten-1)/100;
    nz = get_block([nstem '.fmapqz'],'byte',nst-1,nen-1)/100;
  end
  if P.add_avgnmm
    tm = get_block([tstem '.avgnmm'],'byte',tst-1,ten-1)/10;
    nm = get_block([nstem '.avgnmm'],'byte',nst-1,nen-1)/10;
  end
  if P.add_nuwp
    tw = get_block([tstem '.nuwp'],'byte',tst-1,ten-1);
    nw = get_block([nstem '.nuwp'],'byte',nst-1,nen-1);
  end
  fprintf('  end2: assigning\n');
  tpos1 = 1;
  npos1 = 1;
  for j=1:length(idx2), i=idx2(j);
    tpos2 = tpos1+ten(j)-tst(j);
    npos2 = npos1+nen(j)-nst(j);
    if st<=tlen
      if P.add_coverage, X.coverageT2(i) = max(tc(tpos1:tpos2)); end
      if P.add_fmapqz, X.fmapqzT2(i) = max(tz(tpos1:tpos2)); end
      if P.add_avgnmm, X.avgnmmT2(i) = max(tm(tpos1:tpos2)); end
      if P.add_nuwp, X.nuwpT2(i) = max(tw(tpos1:tpos2)); end
    end
    if st<=nlen
      if P.add_coverage, X.coverageN2(i) = max(nc(npos1:npos2)); end
      if P.add_fmapqz, X.fmapqzN2(i) = max(nz(npos1:npos2)); end
      if P.add_avgnmm, X.avgnmmN2(i) = max(nm(npos1:npos2)); end
      if P.add_nuwp, X.nuwpN2(i) = max(nw(npos1:npos2)); end
    end
    tpos1=tpos2+1;
    npos1=npos2+1;
  end
 end  % if ~P.skip_end2

end

% convert negative coverage values (these are on a log scale)

if P.add_coverage
  y = -1:-1:-127;
  x = round(112+exp(2-y/5));
  X.coverageN1(X.coverageN1<0) = x(-X.coverageN1(X.coverageN1<0));
  X.coverageT1(X.coverageT1<0) = x(-X.coverageT1(X.coverageT1<0));
  if ~P.skip_end2
    X.coverageN2(X.coverageN2<0) = x(-X.coverageN2(X.coverageN2<0));
    X.coverageT2(X.coverageT2<0) = x(-X.coverageT2(X.coverageT2<0));
  end
end

