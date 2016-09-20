function p=single_screen_prob(x,n,r,prob_type,L)

if ~exist('prob_type','var')
  prob_type='poisson';
end

p=zeros(size(x));
switch prob_type
 case 'poisson'
  p=poisspdf(x,n*r);
 case 'binom'
  p=binopdf(x,n*L,r/L);
 case 'gaussian'
  if x(1)~=0
    p=normcdf([x(1)-0.5 x+0.5],n*r,sqrt(n*r));
    p=diff(p);
  else
    p=normcdf([-Inf x+0.5],n*r,sqrt(n*r));
    p=diff(p);    
  end
end
