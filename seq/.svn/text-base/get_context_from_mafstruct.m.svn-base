function n = get_context_from_mafstruct(M,ncat,context_field_name)
% function n = get_context_from_mafstruct(M,ncat)
%
% Mike Lawrence 2010-01-27

if ~exist('context_field_name','var'), context_field_name='context'; end
if context_field_name(1)=='/', context_field_name = regexprep(context_field_name,'^.*/([^/]*)$','$1'); end

if isempty(ncat)
  % infer from categs.txt
  tmp = load_struct(['/xchip/cga1/lawrence/db/' context_field_name '/categs.txt'],'%f%s');
  if tmp.num(1)==0, fprintf('Warning: category 0 encountered: was assuming 1-based\n'); end
  ncat = tmp.num(end);
end

if isfield(M,'classification')
  M = reorder_struct(M,grep('SNP',M.classification,1));
end

C = getfield(M,context_field_name);
if ~isnumeric(C), C = str2double(C); end

n = nan(ncat,4);
base='ACGT';

if ~isfield(M,'newbase')
  M.newbase = find_newbase(M);
end

for i=1:4
  q = C(strcmp(M.newbase,base(i)));
  if any(q==0), fprintf('Warning: category 0 encountered: was assuming 1-based\n'); end
  n(:,i) = histc(q,1:ncat);
end
