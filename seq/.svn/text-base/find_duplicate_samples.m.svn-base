function [M ni fi ov Q] = find_duplicate_samples(M,P)
% M = find_duplicate_samples(M,P)
%
% finds samples that are likely to be duplicates, on the basis of shared mutations
%
% Mike Lawrence 2012-10-02

if nargout==4, error('return arguments are different now'); end

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'max_overlap_muts_allowed',10);
P = impose_default_value(P,'max_pct_overlap_allowed',20);

P = impose_default_value(P,'real_recurrent_site_genes_to_exclude',...
   {'NRAS','KRAS','HRAS','TP53','PIK3CA','SPOP','SF3B1','U2AF1','GATA3','AKT1',...
   'BRAF','APC','ATM','FBXW7','SMAD4','SMAD2','ERBB3','MYD88','IDH1','IDH2','IDH3',...
   'EGFR','NFE2L2','KEAP1','CDKN2A','DNMT3A','FLT3','KIT','VHL','PBRM','CTNNB1','ARID1A',...
   'PTEN','NF1','RB1'});

if isfield(M,'Reference_Allele') || isfield(M,'ref_allele')
  tmp=M;
  M=[];
  M.mut = tmp;
  clear tmp;
end

if isfield(M,'patient') && ~isfield(M,'pat'), M=rename_field(M,'patient','pat'); end

if ~isfield(M,'pat')
  [M.pat.name ui M.mut.pat_idx] = unique(M.mut.patient);
end

demand_fields(M,{'pat','mut'});
np = slength(M.pat);

if isfield(M.mut,'gene')
  is_excluded = ismember(M.mut.gene,P.real_recurrent_site_genes_to_exclude);
  fprintf('Excluding %d mutations in driver genes\n',sum(is_excluded));
elseif isfield(M.mut,'Hugo_Symbol')
  is_excluded = ismember(M.mut.Hugo_Symbol,P.real_recurrent_site_genes_to_exclude);
  fprintf('Excluding %d mutations in driver genes\n',sum(is_excluded));
else
  is_excluded = false(slength(M.mut),1);
  fprintf('NOTE: no gene information found!  Cannot exclude driver genes\n');
end

demand_fields(M.mut,{'chr','start'});
if ~isnumeric(M.mut.chr), M.mut.chr=convert_chr(M.mut.chr); M.mut = reorder_struct(M.mut,~isnan(M.mut.chr)); end
M.mut = make_numeric(M.mut,'start');

if isfield(M.mut,'pat_idx')
  % great
elseif isfield(M.mut,'patient') && isfield(M,'pat') && isfield(M.pat,'name')
  M.mut.pat_idx = listmap(M.mut.patient,M.pat.name);
else
  error('missing per-mutation patient information');
end
if isfield(M.mut,'newbase_idx')
  % great
elseif isfield(M.mut,'newbase')
  M.mut.newbase_idx = listmap(M.mut.newbase,{'A','C','G','T'});
else
  error('missing newbase information');
end

w = cell(np,1);
widx = cell(np,1);
for p=1:np
  idx = find(M.mut.pat_idx==p);
  w{p} = double(1e13*M.mut.chr(idx)+100*M.mut.start(idx)+M.mut.newbase_idx(idx));
  w{p}(is_excluded(idx)) = nan;  % will fail to match in ismember()
  widx{p} = idx;
end

M.pat.dup_partner = repmat({[]},np,1);
M.pat.dup_nmuts = repmat({[]},np,1);
M.pat.dup_fracmuts = repmat({[]},np,1);

M.mut.dup_count = zeros(slength(M.mut),1);

ni = nan(np,np); fi = nan(np,np); ov=false(np,np);
for p1=1:np-1, if ~mod(p1,10), fprintf('%d/%d ',p1,np); end
  for p2=p1+1:np
    n1 = length(w{p1});
    n2 = length(w{p2});
    midx1 = widx{p1}(ismember(w{p1},w{p2}));
    midx2 = widx{p2}(ismember(w{p2},w{p1}));
    midx = union(midx1,midx2);
    M.mut.dup_count(midx) = M.mut.dup_count(midx) + 1;
    ni(p1,p2) = length(midx);
    fi(p1,p2) = ni(p1,p2)/(n1+n2-ni(p1,p2));
    ov(p1,p2) = (ni(p1,p2)>=P.max_overlap_muts_allowed || fi(p1,p2)>=(P.max_pct_overlap_allowed/100));
    if ov(p1,p2)
      fprintf('\n(%d)  %s   %d muts     (%d)  %s   %d muts    overlap = %d muts (%.2f%%)\n',...
              p1,M.pat.name{p1},n1,p2,M.pat.name{p2},n2,ni(p1,p2),fi(p1,p2)*100);
      M.pat.dup_partner{p1}(end+1) = p2;
      M.pat.dup_partner{p2}(end+1) = p1;
      M.pat.dup_nmuts{p1}(end+1) = ni(p1,p2);
      M.pat.dup_nmuts{p2}(end+1) = ni(p1,p2);
      M.pat.dup_fracmuts{p1}(end+1) = fi(p1,p2);
      M.pat.dup_fracmuts{p2}(end+1) = fi(p1,p2);
    end
  end
end, fprintf('\n');

M.pat.ndup_partners = cellfun(@length,M.pat.dup_partner);

% print list of sample "cliques"
try
  Q={};used = false(np,1);
  for i=1:np, if used(i), continue; end
    clique = i; done = false;
    while(~done); done = true;
      for j=1:length(clique), k = clique(j);      
        add = setdiff(find((ov(k,:)'|ov(:,k)) & ~used),clique);
        if ~isempty(add), done = false; clique = union(clique,add); end
      end
    end
    if length(clique)>1, Q{end+1} = as_column(M.pat.name(clique)); end
    used(clique) = true;
  end
  for i=1:length(Q)
    fprintf('Clique %d (%d members)\n',i,length(Q{i}));
    pr(M.pat,Q{i});
    fprintf('\n');
  end
catch me
  fprintf('Problem in clique finding...\n');
  keyboard
end



