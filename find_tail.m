function t=find_tail(val,dat,tail)
% t=find_tail(vals,dat,tail)
%    returns in t a vector of the same size as vals
%    where in each element are the number of enties in dat that 
%    are >= (tail==1) or <= (tail==-1) in dat
%    val and dat are assumed to be column vectors
% 
% Gaddy Getz
% Cancer Genomics
% The Broad Institute
% gadgetz@broad.mit.edu
%


ndat=size(dat,1);
nv=size(val,1);
sdat=sort(dat);
%val=val-tail*eps;
[sval,svali]=sort(val);
[tmp,svalrevi]=sort(svali);
if tail==-1 
  v=[ ones(ndat,1); zeros(nv,1)];
  [tmp,ord]=sort([ dat; sval ]);
  t=cumsum(v(ord));
elseif tail==1
  v=[ zeros(nv,1); ones(ndat,1)];
  [tmp,ord]=sort([ sval; dat ]);
  t=ndat-cumsum(v(ord));
else
  error('tail should be -1 or 1');
end
t=t(find(v(ord)==0));
t=t(svalrevi);

