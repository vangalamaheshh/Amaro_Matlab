p = 10;
n = 10000;
 
rls = recursiveLS(p);
 
x = randn(n,p);
b = 3*randn(p,1);
y = x*b + 0.1*randn(n,1);
 
thetas = zeros(length(x),p);
for i = 1 : length(x)
    
    thetas(i,:) = step(rls, y(i), x(i,:));
end