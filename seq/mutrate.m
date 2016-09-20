function [ratio sd ntot Ntot] = mutrate(C,x1,y1,z1,x2,y2,z2)

if nargin~=4 && nargin~=7, error('takes 4 or 7 input arguments'); end

ntot = fullsum(C.n.gcn(x1,y1,z1));
Ntot = fullsum(C.cov.gc(x1,y1));
if nargin==7
  ntot = ntot + fullsum(C.n.gcn(x2,y2,z2));
  Ntot = Ntot + fullsum(C.cov.gc(x2,y2));
end
[ratio sd] = ratio_and_sd(ntot,Ntot);
fprintf('%0.2f +- %0.2f\n',ratio*1e6,sd*1e6);

if nargout==0
  clear ratio sd ntot Ntot
end
