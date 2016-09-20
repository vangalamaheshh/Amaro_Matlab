function D=make_D(dat,gacc,gdesc,sdesc)

D.dat=dat;
if ~exist('gacc','var') || isempty(gacc)
  gacc=cellstr(num2str((1:size(dat,1))','%d'));
end
D.gacc=gacc;
if  ~exist('gdesc','var') ||isempty(gdesc)
  D.gdesc=gacc;
else
  D.gdesc=gdesc;
end
if ~exist('sdesc','var') || isempty(sdesc)
  sdesc=cellstr(num2str((1:size(dat,2))','%d'));
end

D.sdesc=sdesc;
