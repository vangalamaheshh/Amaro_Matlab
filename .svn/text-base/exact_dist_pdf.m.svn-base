function d=exact_dist_pdf(n,f,range,eps,maxn)

if ~exist('eps','var')
  eps=1e-100;
end

m=size(n,2); % m types
v=zeros(length(range),m);

for i=1:m
%  x1=[1 1-binocdf(0:(n(1,i)-1),n(1,i),f(i))];
%  y1=-log(x1);
  if exist('maxn','var')
    x=ln_binopdf(0:maxn,n(1,i),f(i));
  else
    x=ln_binopdf(0:n(1,:),n(1,i),f(i));
  end
  y=-x; % ln_binocdf(0:n(1,i),n(1,i),f(i),1);

  large=find(x>log(eps));
  for j=1:length(large)
%    bin=floor((y(large(j))-range(1))/(range(2)-range(1)))+1;
%    if bin<1 || bin>length(range)
%      bin=[];
%    end
    bin=find(range(1:(end-1))<=y(large(j)) & range(2:end)>y(large(j)));
%     disp([ y(large(j)) exp(-y(large(j))) bin]);

%    disp([ length(bin) length(bin1) bin bin1]);
    if ~isempty(bin)
      v(bin,i)=v(bin,i)+exp(x(large(j)));
%      disp([ bin y(large(j)) v(bin,i) exp(x(large(j)))]);
    end
  end
%   disp(i);
end
if (0)
  dd=v(:,1);
  for i=2:m
    dd=conv(dd,v(:,i));
    dd=dd(1:min(length(dd),2*length(range)));
  end
  d=dd;
else
  d=conv_many_fft(v);
end




