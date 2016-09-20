function s=intersect_many(varargin)

s=[];
if ~isempty(varargin)
  s=varargin{1};
  if strcmp(varargin{end},'rows')
    for i=2:(nargin-1)
      s=intersect(s,varargin{i},'rows');
    end
  else
    for i=2:nargin
      s=intersect(s,varargin{i});
    end
  end    
else
    s=[];
end
