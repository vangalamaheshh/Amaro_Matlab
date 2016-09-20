function a = isodd(b)

if ~isnumeric(b)
  a=nan(size(b));
else
  a=mod(b,2);
end
