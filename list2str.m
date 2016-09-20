function str=list2str(tbl,sep,empty_element)

if ~exist('sep','var') || isempty(sep)
  sep={char(9),sprintf(newline)};
end
if ~exist('empty_element','var')
  empty_element='---';
end

c=1;
r=1;
str=[];
for r=1:size(tbl,1)
  for c=1:size(tbl,2)
    x=tbl{r,c};
    if isempty(x)
      str=[ str sprintf('%s',empty_element)];
    elseif ischar(x)
      str=[ str sprintf('%s',x)];
    elseif isnumeric(x)
      if prod(size(x))==1
        if x==round(x)
          str=[str sprintf('%d',x)];
        else
          str=[str sprintf('%f',x)];
        end
      else
        x=mat2cell(x,ones(1,size(x,1)),ones(1,size(x,2)));
        str=[str list2str(x,sep(2:end,:))];
      end
    elseif iscell(x)
      str=[str list2str(x,sep(2:end,:))];
    else
      error('unexpected type');
    end
    if c<size(tbl,2)
      str=[str sprintf('%s',sep{1,1})];
    end
  end
  if c<size(tbl,1)
    str=[str sprintf('%s',sep{1,2})];
  end
end
