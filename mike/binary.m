function s = binary(d,width)
%
% binary(d)
%
% returns string representation of binary expansion of decimal integer d
%
% width (optional) = number of bits (default is lowest possible multiple of 8)
%
% Mike Lawrence 2008-08-21
%

if ischar(d), d=double(d); end

if length(d)==1
  if ~exist('width','var'), width = 8*ceil(log2(d+1)/8); end
  s = get_binary(d,width);
else
  if ~exist('width','var'), width = 8*ceil(log2(fullmax(d)+1)/8); end
  s = cell(size(d));
  for i1=1:size(d,1)
  for i2=1:size(d,2) 
  for i3=1:size(d,3)
  for i4=1:size(d,4)
  for i5=1:size(d,5)
    s{i1,i2,i3,i4,i5} = get_binary(d(i1,i2,i3,i4,i5),width);
  end,end,end,end,end
end

  function s = get_binary(d,width)
    b = (bitand(d,2.^((width-1):-1:0))>0);
    s = num2str(b);
    s(end-23:-24:1) = '.';
    s(s==' ')=[];
    s(s=='.')=' ';
  end

end % main function
