function s=params_struct(varargin)

s=[];
for i=1:2:length(varargin)
  if ~ischar(varargin{i})
    disp(['Expected textual argument in arg no. ' num2str(i)]);
  else
    s=setfield(s,varargin{i},varargin{i+1});
  end
end
