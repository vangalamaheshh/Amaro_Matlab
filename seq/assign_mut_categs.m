function [categ categ_ignoring_null_categ] = assign_mut_categs(M,K,P)
% categ = assign_mut_categs(M,K)
%
% Given a set of k categories as struct K with the following fields:
%   left   = subset of 'ACGT', representing 5' base
%   right  = subset of 'ACGT', representing 3' base
%   from   = subset of 'AC', representing mutated base (after strand collapse)
%   change = subset of 'tfs', representing Transition, Flip transversion, Skew transversion
%   name   = "indel", "null", "indel+null", "double_null" or name of point-mutation category
%
% and a set of mutations M with
%   classification = SNP, DNP, Indel/Ins/Del, or Complex_substitution
%         SNP, DNP = added to the "point" category according to "context65"
%         Indel, Complex_substitution = added to "indel" category
%
%   context65 = 65-category context from get_context
%   newbase (or ref_allele, tum_allele1, and tum_allele2)
%
%   (if "null" assignments requested), also need: type
%   (if "double_null" assignments requested), also need: gene, patient
%
% Returns categ, telling which of the k categories the mutation should be counted in
%   (zero indicates a mutation is not to be counted toward any category)
%
% Mike Lawrence 2010-01-27

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'add_double_null_category',false);

fprintf('Assigning mutation categories...\n');

if isfield(M,'mut'), error('first argument should be M.mut, not M'); end

if ~isfield(M,'context65') && isfield(M,'context')
  fprintf('No "context65" field found in M: will try using "context" field\n');
  M.context65 = M.context;
end

require_fields(M,{'context65','classification','type','newbase'});

require_fields(K,{'left','right','from','change','name'});
nk = slength(K);

M = make_numeric(M,'context65');
if mean(M.context65>=1 & M.context65<=65)<0.9, error('M.context65 is incomplete or damaged'); end
%if max(M.context65)<10, error('M.context65 needs to be 65-category context'); end

% fix blank newbase
if any(cellfun('isempty',M.newbase))
  M.newbase = find_newbase(M);
end

% fix blank classification
idx = find(cellfun('isempty',M.classification));
if ~isempty(idx)
  fprintf('%d/%d entries have blank M.classification: assuming "SNP" for all of these.\n',length(idx),slength(M));
  M.classification(idx) = repmat({'SNP'},length(idx),1);
end

idx1 = grepi('^(SNP|DNP|TNP|ONP)$',M.classification,1);   % "non-indel"
idx2 = grepi('^(Indel|Ins|Del|Consolidated|Complex_substitution)$',M.classification,1);   % "indel/complex"
idx3 = setdiff((1:slength(M)),[idx1;idx2]);
if ~isempty(idx3)  
  fprintf('Unknown M.classification:\n');
  disp(unique(M.classification(idx3)));
%  error('Please correct and retry.');
  fprintf('Will classify these as Indel/Complex.  If this is not correct, please update assign_mut_categs.m\n');
  idx2 = union(idx2,idx3);
end

double_null_categ = grepi('null2|double.?null',K.name,1);
if length(double_null_categ)>1, error('Only one "double_null" category allowed'); end

null_categ = grepi('null',K.name,1);
null_categ = setdiff(null_categ,double_null_categ);
if length(null_categ)>1, error('Only one "null" category allowed'); end

indel_categ = grepi('indel',K.name,1);
indel_categ = setdiff(indel_categ,double_null_categ);
if length(indel_categ)>1, error('More than one indel category not supported'); end
if length(indel_categ)<1 && length(idx2)>0
  fprintf('WARNING:  No indel category specified; will ignore indels in the MAF!');
  ignore_indels = true;
else
  ignore_indels = false;
end

c = assign_65x4_to_categ_set(K);
c = bsxfun(@times,c,1:slength(K));
if length(indel_categ>0), c(c==indel_categ)=0; end  % indel categories not in consideration here
if length(null_categ>0), c(c==null_categ)=0; end  % null categories not in consideration here
if length(double_null_categ>0), c(c==double_null_categ)=0; end  % null categories not in consideration here
c = squeeze(max(c,[],2));

base = {'A','C','G','T'};
newbase = listmap(regexprep(M.newbase,'^(.).*','$1'),base);

nm = slength(M);
categ = zeros(nm,1);
for newbasei=1:4
  idx = idx1(newbase(idx1)==newbasei & ~isnan(M.context65(idx1)));
  categ(idx) = c(M.context65(idx),newbasei);
end
if ignore_indels
  categ(idx2) = 0; % not in any category
else
  categ(idx2) = indel_categ;
end
num_weird_cases = sum(isnan(M.context65(idx1)) | isnan(newbase(idx1)));

% OLDER SLOW LOOP METHOD
%num_weird_cases = 0;
%ismember1 = ismember(1:nm,idx1);   % "non-indel"
%ismember2 = ismember(1:nm,idx2);   % "indel/complex"
%for i=1:nm
%  if ismember1(i)
%    if isnan(M.context65(i)) || isnan(newbase(i))
%      % weird case
%      num_weird_cases = num_weird_cases + 1;
%      categ(i) = 0;  % category undetermined
%    else
%      categ(i) = c(M.context65(i),newbase(i));
%    end
%  elseif ismember2(i) % indel
%    if ignore_indels
%      categ(i) = 0;  % not in any category
%    else
%      categ(i) = indel_categ;
%    end
%  else error('Inconsistent behavior in assign_mut_categs.m');
%  end
%end

if num_weird_cases/slength(M) >= 0.1 || (num_weird_cases/slength(M) >= 0.01 && num_weird_cases>100)
  fprintf('*************************************************************************************\n');
  fprintf('WARNING:   POSSIBLE BUILD MISMATCH BETWEEN context65_dir AND MUTATION ANNOTATIONS    \n');
  fprintf('*************************************************************************************\n');
  keyboard
end

categ_ignoring_null_categ = categ;

% handle "null" and "double-null" categories if provided for
if ~isempty(null_categ) || ~isempty(double_null_categ)
  flds = fieldnames(M);
  typefield = grepi('^(type|Variant_Classification)$',flds);
  if isempty(typefield)
    error('adding null category: can''t find type information');
  end
  Mtype = getfield(M,typefield{1});

  null_idx = grepi(['Frame.?shift|Non.?sense|Read.?through|Non.?stop|'...
                  'Splice|Stop.?codon|Init.?Met|De.?Novo.?Start'],Mtype,1);
  if ~isempty(null_categ)
    categ(null_idx) = null_categ;
  end

  if ~isempty(double_null_categ)
    patfield = grepi('^(patient|patient_name|Tumor_Sample_Barcode)$',flds);
    genefield = grepi('^(gene|gene_name|Hugo_Symbol)$',flds);
    if isempty(patfield) || isempty(genefield)
      error('double_null: can''t find patient+gene information');
    end
    Mpat = getfield(M,patfield{1});
    Mgene = getfield(M,genefield{1});

    [tmp tmp pj] = unique(Mpat);
    [tmp tmp gj] = unique(Mgene);

    is_coding = nan(slength(M),1);
    is_coding(grepi('missense|nonsense|silent|splice|synon|frame|shift|non.?stop|read.?thr',Mtype,1)) = 1;
    is_coding(grepi('intron|utr|promoter|flank',Mtype,1)) = 0;
    is_coding(grepi('igr',Mtype,1)) = nan;    % nan's will all segregate individually --> no IGR double_null's
    is_coding(strcmpi('',Mgene)|strcmpi('None',Mgene)|strcmpi('Unknown',Mgene)) = nan;

    [u tmp uj] = unique([pj(null_idx) gj(null_idx) is_coding(null_idx)],'rows');
    h = histc(uj,1:length(u));
    dn = find(h>=2);
    n_dn = 0; n_ua = 0;
    for i=1:length(dn)
      idx = null_idx(uj==dn(i));
      categ(idx(1)) = double_null_categ;    % mark the first as "double_null"
      categ(idx(2:end)) = 0;                % mark the other(s) as "unassigned"
      n_dn = n_dn + 1;
      n_ua = n_ua + length(idx)-1;
    end
    fprintf('Designated %d double-null mutations.\n',n_dn);
    fprintf('Removed assignments of their %d partner-null mutations; these will be omitted from analysis.\n',n_ua);
  end
end

% summarize disposition of mutations
categname = nansub(num2cellstr(0:slength(K)),categ+1,'?');
xcount(M.type,categname)
try
  if any(categ==0)
    pr(0:slength(K),['unassigned';K.name]);
  else
    pr(1:slength(K),K.name);
  end
catch me
  fprintf('Problem printing list of categories\n');
end

% final quality check
weird_cases = find(categ==0 | categ_ignoring_null_categ==0);
if ignore_indels, weird_cases = setdiff(weird_cases,idx2); end
num_weird2_cases = length(weird_cases);
if num_weird2_cases/slength(M) >= 0.1 || (num_weird2_cases/slength(M) >= 0.01 && num_weird2_cases>100)
  fprintf('**************************************************************************************\n');
  fprintf('WARNING:   POSSIBLE BUILD MISMATCH.  Too many collisions between context65 and newbase\n');
  fprintf('**************************************************************************************\n');
  keyboard
end
