function C = load_coverage(C,P)
% Mike Lawrence 2009-2010

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'col1_is_patient_name',true);
P = impose_default_value(P,'impute_full_coverage',false');

require_fields(C,{'ns','ncat','sample','file'});

% make sure all coverage files exist
if ~P.impute_full_coverage
  covfile = cell(C.ns,1);
  for i=1:C.ns
    if isfield(C.sample,'covfile')
      covfile{i} = C.sample.covfile{i};
    else
      covfile{i} = ['/cga/tcga-gsc/home/lawrence/' C.sample.dir{i} '/' C.file.cov];
    end
    demand_file(covfile{i});
  end
end

% load target file
C.targ = load_target_file(C.file.targ);
fprintf('Sorting target list.\n');
C.targ = sort_struct(C.targ,{'chr','start','end'});
C.nt = slength(C.targ);

C.cov = nan(C.nt,C.ns,C.ncat);
if ~P.impute_full_coverage
  % load coverage
  for i=1:C.ns
    fprintf('%d/%d Loading %s\n',i,C.ns,covfile{i});
    
    fmt = '%s';
    if P.col1_is_patient_name, fmt = [fmt '%s']; end
    fmt = [fmt repmat('%f',1,C.ncat+3)];
    tbl = read_table(covfile{i},fmt,char(9),0);
    if size(tbl.dat{1},1)~=size(C.cov,1), fprintf('\tWrong number of exons!\n',i,covfile{i}); keyboard; end
    
    X = [];
    X.chr = tbl.dat{2+P.col1_is_patient_name};
    X.start = tbl.dat{3+P.col1_is_patient_name};  
    X.end = tbl.dat{4+P.col1_is_patient_name};
    [X ord] = sort_struct(X,{'chr','start','end'});
    if any(X.chr~=C.targ.chr | X.start~=C.targ.start | X.end~=C.targ.end), fprintf('\tWrong exon coordinates!\n'); keyboard; end
    for c=1:C.ncat, C.cov(:,i,c) = tbl.dat{4+P.col1_is_patient_name+c}(ord); end
  end
else   % P.impute_full_coverage
  fprintf('Imputing full coverage...\n');
  for chr=1:24,fprintf('chr%d: ',chr);
    k = load([C.file.categdir '/chr' num2str(chr) '.mat']);
    flds = fieldnames(k);
    k = getfield(k,flds{1});
    idx = find(C.targ.chr==chr); nidx=length(idx);
    for ti=1:nidx
      if ~mod(ti,1000), fprintf(' %d/%d',ti,nidx); end
      i=idx(ti);
      C.cov(i,1,:) = histc(k(C.targ.start(i):C.targ.end(i)),1:C.ncat);
    end, fprintf('\n');
  end
  for s=2:C.ns, C.cov(:,s,:) = C.cov(:,1,:); end
end

