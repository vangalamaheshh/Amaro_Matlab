function s=num2chromosome(n)
% NUM2CHROMOSOME convert chromosome index to string representation
%   S = NUM2CHROMOSOME(N)
%   Returns a cell array of strings corresponding the chromosome
%   numbers in N. If N is a scalar, a simple string is returned.
%

% GISTIC software version 2.0
% Copyright (c) 2011 Gad Getz, Rameen Beroukhim, Craig Mermel, 
% Jen Dobson, Steve Schumacher, Nico Stransky, Mike Lawrence, 
% Gordon Saksena, Michael O'Kelly, Barbara Tabak
% All Rights Reserved.
%
% See the accompanying file LICENSE.txt for licensing details.

conv={'1','2','3','4','5','6','7','8','9','10',...
      '11','12','13','14','15','16','17','18','19','20',...
      '21','22','X','Y','M','XY'};
% undo SegArray
n = full(n);

if numel(n)>1
  s=conv(n);
%   s=cell(size(n));
%   for i=1:size(s,1)
%     for j=1:size(s,2)
%       s{i,j}=num2chromosome(n(i,j));
%     end
%   end
else
  s=conv{n};  
%   if n<=22
%     s=num2str(n);
%   elseif n==23
%     s='X';
%   elseif n==24
%     s='Y';
%   else
%     s='M'
%   end
end

