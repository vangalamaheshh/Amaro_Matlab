function h = hist3d(a,b,c,firsta,lasta,firstb,lastb,firstc,lastc)

if nargin~=9, error('requires 9 input arguments'); end

h = zeros(lasta-firsta+1,lastb-firstb+1,lastc-firstc+1);

for bi=firstb:lastb, for ci=firstc:lastc
  h(:,bi,ci) = histc(a(b==bi & c==ci),firsta:lasta);
end,end

