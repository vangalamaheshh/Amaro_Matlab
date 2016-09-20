function D=threshold(D,val)

D=add_history(D,mfilename,val);

if isfield(D,'tdat')
  warning([ mfilename ': obsolete']);
  D=rmfield(D,'tdat');
end

%D.tdat=D.dat;
%D.tdat(D.dat<val)=val;

D.dat=D.dat;
D.dat(D.dat<val)=val;


