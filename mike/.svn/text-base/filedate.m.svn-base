function sz = filedate(fname)

d = dir(fname);
sz = nan(length(d),1);
for i=1:length(d)
  sz(i) = d(i).date;
end
