function [u,ui,uj]=unique_keepord(x,isrows)


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

if exist('isrows','var') && strcmp(isrows,'rows')
  [uu,ui,uj]=unique(x,'rows','first');
%  for i=1:length(ui)
%    o(i)=min(find(uj==i));
%  end
  o=ui;
  
  [so,si]=sort(o);
  u=uu(si,:);
  ui=ui(si);
  [dum,revsi]=sort(si);
  uj=revsi(uj);
else
  [uu,ui,uj]=unique(x,'first');
%  for i=1:length(uu)
%    o(i)=min(find(x==uu(i)));
%  end
  o=ui;
  
  [so,si]=sort(o);
  u=uu(si);
  ui=ui(si);
  [dum,revsi]=sort(si);
  uj=revsi(uj);
end
