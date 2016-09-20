function S = load_picard_struct(fname)
% Mike Lawrence 2009-09-28

l = load_lines(fname);
l = l(~strncmp(l,'#',1));
l = l(~cellfun('isempty',l));
if length(l)<2, error('Problem parsing %s',fname); end

f = split(l{1},char(9));
nf = length(f);

l = l(2:end);
nl = length(l);

C = cell(nl,nf);
for i=1:nl
  v = split(l{i},char(9));
  if length(v)~=nf, error('Inconsistent number of fields in %s',fname); end
  C(i,:) = v;
end

S = [];
f = genfieldname(f);
for i=1:nf
  S = setfield(S,f{i},C(:,i));
end

  
