function [nsil nnon] = get_nsil_nnon_from_mafstruct(M,ncat,context_field_name)
% [nsil nnon] = get_nsil_nnon_from_mafstruct(M,ncat,context_field_name)

if ~exist('context_field_name','var'), context_field_name='context'; end

if isempty(ncat)
  % infer from categs.txt
  tmp = load_struct(['/xchip/cga1/lawrence/db/' context_field_name '/categs.txt'],'%f%s');
  if tmp.num(1)==0, fprintf('Warning: category 0 encountered: was assuming 1-based\n'); end
  ncat = tmp.num(end);
end

M = reorder_struct(M,grep('SNP',M.classification,1));

C = getfield(M,context_field_name);
if ~isnumeric(C), C = str2double(C); end

idx = grepi('silent|synonymous',M.type,1);
nsil = histc(C(idx),1:ncat);

idx = grepi('missense|nonsense|read-through',M.type,1);
nnon = histc(C(idx),1:ncat);
