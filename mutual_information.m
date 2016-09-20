function mi=mutual_information(varargin)

if nargin==1
  m=varargin{1};
  m=m+eps;
  s1=sum(m,1);
  s2=sum(m,2);
  n=sum(s1);
  
  p1=s1./n;
  p2=s2./n;
  
  p=m./n;
  
  % sum(sum(p))
  r=(repmat(p1,2,1).*repmat(p2,1,2));
  % sum(sum(r))
  I=log2(p./r);
  
  mi=sum(sum(p.*I));
else
  both=varargin{1};
  % add an option to do many at once
  disp('FIXME');
end
