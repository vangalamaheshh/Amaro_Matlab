function check_multiM_agreement(M)

if ~iscell(M), error('M is not a cell!'); end

nm = length(M);
if nm==1, fprintf('M has only one element: no check required\n'); return; end

% check ng
for i=2:nm
  if M{i}.ng ~= M{1}.ng, error('multiM disagreement: different genecounts'); end
end

% check cov.nt
for i=2:nm
  if M{i}.cov.nt ~= M{1}.cov.nt, error('multiM disagremeent: different targetcounts'); end
end

if isfield(M{1}.cov,'gene')
  % check cov.gene
  for i=2:nm
    if slength(M{i}.cov.gene) ~= slength(M{1}.cov.gene), error('multiM disagreement: different genecounts in cov'); end
  end
end

 
