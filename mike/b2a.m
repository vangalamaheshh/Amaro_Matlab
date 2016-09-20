function a = b2a(b)

a = repmat('.',1,4*length(b));

bases = 'ACGT';

ai=1;
for bi=1:length(b)
  y = b(bi);
  for k=1:4
    x = mod(y,4);
    a(ai) = bases(x+1);
    y = bitshift(y,-2);
    ai=ai+1;
  end
end
