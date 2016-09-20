function C = load_coverage_compactly_and_convert_to_mat(C,K,P)
% C = load_coverage_compactly_and_convert_to_mat(C,K,P)
%
% memory-efficient replacement for load_coverage + collapse_coverage_categories + compute_fraction_coverage
%
% 1. loads targ list
%    -- creates C.orig_cov(ns,[original number of uncollapsed categories]) == for category discovery
%    -- creates C.fcov(nt,ns)  == for capture_covplot (entire territory)
%    -- creates C.fcov_coding(nt,ns)  == for capture_covplot (coding part only)
%    -- creates C.gene,C.ng,C.gene.{name,chr,start,end,len,gc,[type]}
%    -- creates C.gene_totcov(ng,ns)   == for building N table
%    -- creates C.gene_cov(ng,ns,ncat) == for building N table
%    -- creates C.gene_sil_cov(ng,ns,ncat) == for analyzing S/NS ratios (included if orig_cat has effect in it)
%    -- creates C.gene_non_cov(ng,ns,ncat) 
%    -- creates C.gene_flank_cov(ng,ns,ncat), which includes all noncoding coverage
%    -- creates C.gene_orig_cov(ng,orig_ncat) == (included optionally)
%    -- creates C.gene_samp_orig_cov(ng,ns,orig_ncat) == (included optionally)
%
%    tabulates territory:
%    -- creates C.orig_terr([original number of uncollapsed categories]) == for category discovery
%    -- creates C.gene_orig_terr(ng,orig_ncat)
%    -- creates C.gene_terr(ng,ncat) == gene territory per final category
%    -- creates C.gene_totterr(ng)   == total territory per gene
%    -- creates C.gene_sil_terr(ng,ncat) == territory broken down by silent/nonsilent
%    -- creates C.gene_non_terr(ng,ncat)
%    -- creates C.gene_flank_terr(ng,ncat), which includes all specified noncoding coverage
%
% 2. for each covfile
%    -- loads textfile (or binary file if available)
%    -- saves as mat (~1/15 the time to load later and ~1/4 the storage space)
%    -- collapses all targets, adds to C.orig_cov
%    -- collapses all categories, adds to C.fcov (will divide by length at end)
%    -- collapses targets to genes
%       -- collapses all categories, adds to C.gene_totcov
%       -- collapses by K, adds to C.gene_cov
%
% Mike Lawrence 2010-09-29

if ~exist('P','var'), P=[]; end
if isfield(P,'include_g_orig_ncat_table')
  fprintf('parameter "include_g_orig_ncat_table" has been renamed "include_gene_orig_ncat_table"\n');
  P = rename_field(P,'include_g_orig_ncat_table','include_gene_orig_ncat_table');
end
P = impose_default_value(P,'include_gene_orig_ncat_table',false);
P = impose_default_value(P,'include_gene_samp_orig_ncat_table',false);
P = impose_default_value(P,'write_mat_files',true);
P = impose_default_value(P,'impute_full_coverage',false);

require_fields(C,{'sample','file'});
if ~isfield(C,'ns'), C.ns = slength(C.sample); end
require_field(C.sample,'name');
if ~P.impute_full_coverage
  require_field(C.sample,'covfile');
end
require_fields(C.file,{'targ','categdir'});

C.file.categs_txt = [C.file.categdir '/categs.txt'];
C.file.categs_fwb = [C.file.categdir '/all.fwb'];
demand_file(C.file.categs_fwb);
C.cat = load_struct(C.file.categs_txt);
if strcmp(C.cat.num{1},'0'), error('Category list includes a zero category!'); end
if ~isfield(C,'ncat')
  C.ncat = slength(C.cat);
else
  if C.ncat ~= slength(C.cat)
    fprintf('WARNING: C.ncat=%d is incorrect: substituting correct value C.ncat=%d\n',C.ncat,slength(C.cat));
    C.ncat = slength(C.cat);
  end
end

binmode = false;
if ~P.impute_full_coverage
  % see whether to use BINARY MODE
  q = length(grep('\.bin$',C.sample.covfile,1));
  if q>=1 && q<C.ns
    error('mixed binary + text files not currently supported');
  elseif q==C.ns
    binmode = true;
    fprintf('all files are *.bin: will read in BINARY MODE\n');
  else
    tmp = regexprep(C.sample.covfile,'(.*)','$1.bin');
    binmode = true;
    for i=1:length(tmp), if ~exist(tmp{i},'file'), binmode = false; break; end, end
    if binmode
      fprintf('Found all *.txt.bin: will read in BINARY MODE\n')
      C.sample.covfile = tmp;
    end
  end

  demand_file(C.sample.covfile);
  if ~binmode && ~isfield(C,'covmat') && P.write_mat_files==true;
    fprintf('Providing default names for covmat to be saved.\n');
    C.sample.covmat = regexprep(C.sample.covfile,'(.*)','$1.mat');
  end
end

% process original categories
try
  map65 = map_categories_to_65(C.file.categs_txt);
  map65mat = sparse(C.ncat,65); for i=1:C.ncat, map65mat(i,map65(i)) = 1; end
catch me
  disp(me); disp(me.message);
  error('ERROR mapping categs_txt to context65 category set');
end
C.cat.refbase = ceil(map65/16);
C.cat.refbase(C.cat.refbase>4) = nan;
tmp = parse(C.cat.name,'(.) in',{'rb'});
tmp.rbi = listmap(tmp.rb,{'A','C','G','T'});
if any(nansum([tmp.rbi -C.cat.refbase],2))
  fprintf('Unexpected error in processing category lists\n');
end

effect_available = false;
try
  map13 = map_categories_to_effect13(C.file.categs_txt);
  [categ13 table13] = get_effect13_categories_list;
  is_coding = (map13~=1);
  effect_available = true;
catch me
  fprintf('No effect information available in category list: will not be able to test NS/S ratios\n');
end

% process and map to final categories
if ~exist('K','var') || isempty(K)
  vanillacatfile = ['/cga/tcga-gsc/home/lawrence/cga/trunk/matlab/seq/' ...
                    'categs_CpGtransit_otherCGtransit_CGtransver_ATmut_IndelAndNull.txt'];
  fprintf('Using default K: using default:\n%s\n',vanillacatfile);
  K = load_struct(vanillacatfile);
end
require_fields(K,{'left','right','from','change','type'});
map65x4 = assign_65x4_to_categ_set(K);
kmat65 = max(map65x4,[],3);
kmat = map65mat * kmat65;   % (map from orig_cat to final categories via 65)

if effect_available
  % plan how to compute expected NS/S ratio for each final category based on coverage totals
  isnon = nan(C.ncat,slength(K),4);
  for i=1:slength(K), isnon(:,i,:) = table13(map13,:); end
  isref = false(C.ncat,slength(K),4);
  for i=1:4, isref(C.cat.refbase==i,:,i) = true; end
  tmp = map65x4(map65,:,:);
  tmp(isref) = 0;
  znon = tmp; zsil = tmp; zflank = tmp;
  znon(isnon~=1) = 0;
  zsil(isnon~=0) = 0;
  zflank(~isnan(isnon)) = 0;
  znon = sum(znon,3)/3;
  zsil = sum(zsil,3)/3;
  zflank = sum(zflank,3)/3;
end

% load target file
C.targ = load_target_file(C.file.targ);
if ~binmode;
  fprintf('Sorting target list.\n');
  C.targ = sort_struct(C.targ,{'chr','start','end'});
end
C.nt = slength(C.targ);

C.orig_cat = C.cat; C = rmfield(C,'cat');
C.orig_ncat = C.ncat; C = rmfield(C,'ncat');
C.categ = K;
C.ncat = slength(K);

% collapse targets to genes
[C.gene.name ui C.targ.gidx] = unique(C.targ.gene);
C.ng = length(C.gene.name);
one_target_per_gene = (C.nt==C.ng);
if C.nt>1e5 && (abs(C.nt-C.ng)/C.nt)<0.001
  fprintf('Approximately one target per gene (%d targets, %d genes): omitting gene collapse\n',C.nt,C.ng);
  one_target_per_gene = true;
end
if one_target_per_gene
  C.gene = C.targ;
  C.ng = C.nt;
  C.gene = rename_field(C.gene,'gene','name');
  C.targ.gidx = (1:C.nt)';
else
  fprintf('Collapsing targets to genes...\n');
  C.gene.chr = nan(C.ng,1);
  C.gene.start = nan(C.ng,1);
  C.gene.end = nan(C.ng,1);
  C.gene.len = nan(C.ng,1);
  if isfield(C.targ,'gc'), C.gene.gc = nan(C.ng,1); end
  if isfield(C.targ,'type'), C.gene.type = nan(C.ng,1); end
  for g=1:C.ng, if ~mod(g,20000), fprintf('%d/%d ',g,C.ng); end
    idx = find(C.targ.gidx==g);
    C.gene.chr(g,1) = C.targ.chr(idx(1));
    C.gene.start(g,1) = min(C.targ.start(idx));
    C.gene.end(g,1) = max(C.targ.end(idx));
    C.gene.len(g,1) = sum(C.targ.len(idx));
    if isfield(C.targ,'gc'), C.gene.gc(g,1) = weighted_mean(C.targ.gc(idx),C.targ.len(idx)); end
    if isfield(C.targ,'type'), C.gene.type(g,1) = C.targ.type(idx(1)); end
  end, if g>=20000, fprintf('\n'); end
  gmat = sparse(C.nt,C.ng); for i=1:C.nt, gmat(i,C.targ.gidx(i)) = 1; end
end

%%%%%%%%%%%%%%%%%%%%% tabulate territory

fwb = org.broadinstitute.cga.tools.seq.FixedWidthBinary(C.file.categs_fwb);
c = nan(C.nt,C.orig_ncat);
for ti=1:C.nt, if ~mod(ti,10000), fprintf('%d/%d ',ti,C.nt); end
  k = fwb.get(C.targ.chr(ti),C.targ.start(ti),C.targ.end(ti));
  c(ti,:) = histc(k,1:C.orig_ncat);
end, fprintf('\n');
fwb.close();

if effect_available
  % for each target, calculate how much of it is coding territory
  C.targ.len_coding = sum(c(:,is_coding),2);
end

if one_target_per_gene
  g = c;
else % collapse to genes
  g = gmat'*c;
end

C.orig_terr = sum(c,1)';
C.gene_orig_terr = g;
C.gene_totterr = sum(g,2);
C.gene_terr = g*kmat;

if effect_available
  C.gene_sil_terr = g*zsil;
  C.gene_non_terr = g*znon;
  C.gene_flank_terr = g*zflank;

  %% check for possible build mismatch
  tt = fullsum(C.gene_terr);
  tc = fullsum(C.gene_sil_terr+C.gene_non_terr);
  if tt > 50e6
    fprintf('Note: total territory is %d: appears to include noncoding areas\n',tt);
  else
    if tc/tt < 0.50
      fprintf('WARNING:  POSSIBLE BUILD MISMATCH\n');
      fprintf('          Most territory appears to be noncoding in categdir\n');
      keyboard;
    end
  end
end


%%%%%%%%%%%%% tabulate per-sample data

if ~P.impute_full_coverage

  % allocate space
  C.orig_cov = zeros(C.ns,C.orig_ncat);
  C.fcov = zeros(C.nt,C.ns);
  if effect_available
    C.fcov_coding = zeros(C.nt,C.ns);
  end
  C.gene_totcov = zeros(C.ng,C.ns);
  C.gene_cov = zeros(C.ng,C.ns,C.ncat);
  if P.include_gene_orig_ncat_table
    C.gene_orig_cov = zeros(C.ng,C.orig_ncat);
  end
  if P.include_gene_samp_orig_ncat_table
    C.gene_samp_orig_cov = zeros(C.ng,C.ns,C.orig_ncat);
  end
  if effect_available
    C.gene_sil_cov = nan(C.ng,C.ns,C.ncat);
    C.gene_non_cov = nan(C.ng,C.ns,C.ncat);
    C.gene_flank_cov = nan(C.ng,C.ns,C.ncat);
  end

  % process each sample
  for i=1:C.ns
    fprintf('%d/%d',i,C.ns);
    
    fname = C.sample.covfile{i};
    fprintf(' Loading %s',fname);

    if ~binmode  % load coverage textfile
      fmt = ['%s%s' repmat('%f',1,C.orig_ncat+3)];
      try
        tbl = read_table(fname,fmt,char(9),0);
      catch me
        error('\tError reading file!');
      end
      if size(tbl.dat,2)~=(5+C.orig_ncat), error('\tWrong number of columns!'); end
      
      % check to make sure it matches the target list
      if size(tbl.dat{1},1)~=C.nt, error('\tWrong number of exons!'); end
      X = []; X.chr = tbl.dat{3}; X.start = tbl.dat{4}; X.end = tbl.dat{5};
      [X ord] = sort_struct(X,{'chr','start','end'});
      if any(X.chr~=C.targ.chr | X.start~=C.targ.start | X.end~=C.targ.end), error('\tWrong exon coordinates!'); end

      % convert to coverage matrix c(nt,orig_ncat)
      c = nan(C.nt,C.orig_ncat);
      for j=1:C.orig_ncat, c(:,j) = tbl.dat{5+j}(ord); end

      if P.write_mat_files==true
        % save as matfile (15x faster to load later, compared to txt file)
        fprintf('\tSaving %s\n',C.sample.covmat{i});
        save(C.sample.covmat{i},'-v7.3','c');
      end

    else   % binary mode: load binary file
      d = dir(fname);
      sz = d.bytes;
      ct = C.nt .* C.orig_ncat;
      if sz==ct, type='byte';
      elseif sz==ct*2, type='short';
      elseif sz==ct*4, type='int';
      else error('size mismatch: %s',fname);
      end
      try
        c = get_block(fname,type,0,ct-1);
      catch me
        error('\tError reading file!');
      end
      c = reshape(c,C.orig_ncat,C.nt)';
    end

    % collapse all targets, add to C.orig_cov
    C.orig_cov(i,:) = sum(c,1);
    % collapse all categories, add to C.fcov (will divide by length at end)
    C.fcov(:,i) = sum(c,2);

    if effect_available
      % do same for C.fcov_coding, but exclude noncoding territory
      C.fcov_coding(:,i) = sum(c(:,is_coding),2);    
    end

    if one_target_per_gene
      g = c;
    else % collapse to genes
      g = gmat'*c;
    end
    if P.include_gene_samp_orig_ncat_table   % save in C.gene_samp_orig_cov
      C.gene_samp_orig_cov(:,i,:) = g;
    end
    if P.include_gene_orig_ncat_table   % add to C.gene_orig_cov
      C.gene_orig_cov = C.gene_orig_cov + g;
    end
    % collapse to K (via 65); add to C.gene_cov
    C.gene_cov(:,i,:) = g*kmat;
    
    if effect_available
      % compute NS/S numerator and denominator for this gene, per (final) category
      for j=1:C.ncat
        C.gene_sil_cov(:,i,j) = g*zsil(:,j);
        C.gene_non_cov(:,i,j) = g*znon(:,j);
        C.gene_flank_cov(:,i,j) = g*zflank(:,j);
      end
    end

    % collapse all categories, add to C.gene_totcov
    C.gene_totcov(:,i) = sum(g,2);
    
    fprintf('\n');
  end % next sample

  % divide by target lengths to finalize C.fcov and C.fcov_coding
  C.fcov = bsxfun(@rdivide,C.fcov,C.targ.len);
  if effect_available
    C.fcov_coding = bsxfun(@rdivide,C.fcov_coding,C.targ.len_coding);
  end    
  
else  % impute full coverage

  fprintf('Will impute full coverage from territory\n');
  return

  C.orig_cov = repmat(C.orig_terr',C.ns,1);
  C.fcov = ones(C.nt,C.ns);
  if P.include_gene_samp_orig_ncat_table
    C.gene_samp_orig_cov = repmat(reshape(C.gene_orig_terr,[C.ng 1 C.orig_ncat]),1,C.ns);
  end
  if P.include_gene_orig_ncat_table
    C.gene_orig_cov = C.gene_orig_terr .* C.ns;
  end
  C.gene_cov = repmat(reshape(C.gene_terr,[C.ng 1 C.ncat]),1,C.ns);
    
  if effect_available
    C.gene_sil_cov = repmat(reshape(C.gene_sil_terr,[C.ng 1 C.ncat]),1,C.ns);
    C.gene_non_cov = repmat(reshape(C.gene_non_terr,[C.ng 1 C.ncat]),1,C.ns);
    C.gene_flank_cov = repmat(reshape(C.gene_flank_terr,[C.ng 1 C.ncat]),1,C.ns);
  end

  C.gene_totcov = repmat(C.gene_totterr,1,C.ns);

end


