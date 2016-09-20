function S = make_numeric(S,varargin)
%
% for a given struct S,
%   performs str2double on all the fields specified
%
% Mike Lawrence 2009-02-24

fields = {};
for i=1:length(varargin)
  fields = [fields varargin{i}];
end

for i=1:length(fields)
  if isfield(S,fields{i})    
    x = getfield(S,fields{i});
    if isnumeric(x) || islogical(x)
      %    fprintf('Warning: %s is already numeric.\n',fields{i});
    else
      x = str2double(x);
      S = setfield(S,fields{i},x);
    end
  else
    fprintf('No such field: %s\n', fields{i});
  end
end
