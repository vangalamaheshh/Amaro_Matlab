function res=conv_many_fft(X)
% convolute columns of X using fft

if nnz(X<0)==0
  all_pos=1;
else
  all_pos=0;
end

n=size(X,1);
m=size(X,2);
X=[X; zeros((m-1)*n-(m-1),m)];
F=zeros(m*n-(m-1),m);
for i=1:m
  F(:,i)=fft(X(:,i),m*n-(m-1));
end

res=real(ifft(prod(F,2)));
if all_pos
  res(res<0)=0;
end

return

n=5000;
m=7;
X=rand(n,m);
X=X./repmat(sum(X,1),n,1);
tic
  res=conv_many_fft(X);
toc

tic
  C=X(:,1);
  for j=2:m
    C=conv(C,X(:,j));
  end
toc
