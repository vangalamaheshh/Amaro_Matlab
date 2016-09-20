function D=log_transform(D)

D=add_history(D,mfilename);

if isfield(D,'tdat')
  warning([ mfilename ': obsolete']);
  if isfield(D,'ltdat')
    D=rmfield(D,'ltdat');
  end
  D.ltdat=log2(D.tdat);
end

D.dat=log2(D.dat);



