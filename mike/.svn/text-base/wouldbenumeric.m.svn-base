function w = wouldbenumeric(x)

if isnumeric(x)
  w = true;
else
  z = mean(~isnan(str2double(x)));
  w = (z>0.5);
end
