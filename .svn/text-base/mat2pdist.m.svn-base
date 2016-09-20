function y = mat2pdist( d )

[n n]=size(d);
y = zeros(1,n*(n-1)/2);
for i=1:n-1
  y( ((i-1)*(2*n-i)/2+1):((i-1)*(2*n-i)/2+n-i) ) = d( i, (i+1):n );
end
