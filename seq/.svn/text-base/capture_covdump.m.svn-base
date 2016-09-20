function capture_covdump(C,P)
% function capture_covdump(C,P)
%
% Mike Lawrence 2009-06-09
%
% now obsolete: happens automatically in capture_covplot (if P.dump_filename provided)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'unit','target');
P=impose_default_value(P,'genelist',[]);
P=impose_default_value(P,'dump_filename','*required*');

fieldname = genfieldname(regexprep(C.sample.med,'-','_'));

if strcmp(P.unit,'target')
  if isempty(P.genelist), idx = 1:C.nt;
  else idx = find(ismember(C.targ.gene,P.genelist)); end
  T = reorder_struct(C.targ,idx);
  for i=1:C.ns, fprintf(' %d/%d',i,C.ns);
    T = setfield(T,fieldname{i},C.cov(idx,i));
  end, fprintf('\n');
elseif strcmp(P.unit,'gene')
  if isempty(P.genelist), idx = 1:C.ng;
  else idx = find(ismember(C.gene.name,P.genelist)); end
  T = reorder_struct(C.gene,idx);
  for i=1:C.ns, fprintf(' %d/%d',i,C.ns);
    T = setfield(T,fieldname{i},C.gcov(idx,i));
  end, fprintf('\n');
else
  error('P.unit must be "target" or "gene"');
end

save_struct(T,P.dump_filename);



