function l=findstrings(strs,s,fld)
% GG Feb 9,2004: works with cell array
% GG Feb 19,2004: if fld is present then looks at field fld of the cell

if ~iscell(strs)
  l=[];
  for i=1:size(strs,1)
    if (findstr(strs(i,:),s))
      l=[l i];
    end
  end
elseif nargin==3
  l=[];
  for i=1:length(strs)
    if (findstr(getfield(strs{i},fld),s))
      l=[l i];
    end
  end
else
  l=[];
  for i=1:length(strs)
    if (findstr(strs{i},s))
      l=[l i];
    end
  end  
end
