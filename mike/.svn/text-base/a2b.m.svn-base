function b=a2b(a)

if mod(length(a),4)~=0, error('input must be in multiples of 4 bytes'), end

a=lower(a);
base = 'acgt';
b=zeros(length(a)/4,1);

bi=1;
ai=1;
j=1;
byte=0;
while ai<=length(a)
  idx = find(base==a(ai));
  if isempty(idx), error('illegal character'); end
  byte=byte/4;
  byte=byte+(4*4*4*(idx-1));
  ai=ai+1;
  j=j+1;
  if j>4
    b(bi)=byte;
    bi=bi+1;
    j=1;
    byte=0;
  end
end

