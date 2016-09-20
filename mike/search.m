function P = search(M,v)



s = size(M);
ndim = length(s);
if ndim>7, error('max dims = 7'), end


sz = ones(7,1);
sz(1:ndim) = s;

P = [];

x=0;
for i1=1:sz(1)
for i2=1:sz(2)
for i3=1:sz(3)
for i4=1:sz(4)
for i5=1:sz(5)
for i6=1:sz(6)
for i7=1:sz(7)
  if M(i1,i2,i3,i4,i5,i6,i7)==v
    coords = [i1 i2 i3 i4 i5 i6 i7];   
    x=x+1;
    P(x,1:length(sz)) = coords(1:length(sz));
  end
end,end,end,end,end,end,end

if x>0
  P = P(:,1:ndim);
else
  P = [];
end
