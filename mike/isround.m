function b = isround(a)

if ~isnumeric(a)
  b = false;
else
  b = (round(a)==a);
end
