function param=enum_param(param,enums)
% param=enum_param(param,enums)
% converts a number or a string of a number to a string accrodign to enums
% enums={val1,'string1';val2,'string2';...}

if ischar(param)
  param_num=str2num(param);
else
  param_num=param;
end
if ~isempty(param_num)
  i=find(cat(1,enums{:,1})==param_num);
  if isempty(i)
    error('No such value allowed');
  end
  param=enums{i,2};
else
  error('param is not a number or a string of a number');
end
