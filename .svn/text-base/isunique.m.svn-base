function res=isunique(x,do_rows)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if size(x,1)==1
  x=x';
end

if exist('do_rows','var')
  [u,ui,uj]=unique(x,'rows');
else
  [u,ui,uj]=unique(x);
end

if size(u,1)==size(x,1)
  res=true;
else
  res=false;
end
