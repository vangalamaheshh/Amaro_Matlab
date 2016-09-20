function [Mt,m1,m2]=match_num_sets(set1,set2);


% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

Mt=[];
[uset1,ui1,uj1]=unique(set1);
[uset2,ui2,uj2]=unique(set2);

% check if both are unique
if length(uset1)==length(set1) && length(uset2)==length(set2) 
  
  c=intersect(set1,set2);
  
  m2=find(ismember(set2,c));
  x=find(ismember(set1,c));
  
  [s1,s1i]=sort(set2(m2));
  [tmp,rs1i]=sort(s1i);
  
  [s2,s2i]=sort(set1(x));
  [tmp,rs2i]=sort(s2i);
  
  m1=x(s2i(rs1i));
  Mt=sparse(m1,m2,true,length(set1),length(set2));
else
  [tmp,um1,um2]=match_num_sets(uset1,uset2);
  uMt=sparse(um1,um2,ones(length(um1),1),length(uset1),length(uset2));
  Mt=uMt(uj1,uj2);
  [m1,m2,tmp]=find(Mt);
end
