function res=convert_enum(str,enums)
if prod(size(enums{1,2}))>1
  use_cell=1;
  res={};
else
  use_cell=0;
  res=NaN;
end
if ~ischar(str) && length(str)>1
  if iscell(str)
    for i=1:length(str)
      if use_cell
        res{i}=convert_one_enum(str{i},enums);
      else
        res(i)=convert_one_enum(str{i},enums);
      end
    end
  else
    for i=1:length(str)
      if use_cell
        res{i}=convert_one_enum(str(i),enums);
      else
        res(i)=convert_one_enum(str(i),enums);
      end
    end
  end
else
  res=convert_one_enum(str,enums);
end

function res=convert_one_enum(str,enums)  
res=NaN;
if ischar(str)
  for i=1:length(enums)
    if strcmp(str,enums{i,1})
      res=enums{i,2};
      break;
    end
  end
else
  for i=1:length(enums)
    if str==enums{i,1} % GG (3/14/06): strcmp(str,enums{i,1})
      res=enums{i,2};
      break;
    end
  end
end

