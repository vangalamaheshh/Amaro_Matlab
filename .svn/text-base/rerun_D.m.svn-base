function D1=rerun_D(D,s)

if ~isfield(D,'history')
  error('D has no history');
end

if ~exist('s','var')
  s=length(D.history);
else
  s=min(s,length(D.history));
end

D1=D;
D1.dat=D.orig;
D1=rmfield(D1,'history');

for i=1:s
  func_name=D.history{i}{2};
  params=D.history{i}(3:end);
  D1=feval(func_name,D1,params{:});
end
