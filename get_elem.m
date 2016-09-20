function y=get_elem(x,varargin)

st='y=x(';
if ~isempty(varargin)
  for i=1:(length(varargin)-1)
    st=[st 'varargin{' num2str(i) '},'];
  end
  st=[st 'varargin{' num2str(length(varargin)) '});'];
  eval(st);
else
  y=x;
end

