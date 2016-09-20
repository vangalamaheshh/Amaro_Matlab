%1) log likelihood poisson regression
% L(lambda|x)= prod(i=1...n)(lambda^x(i)*exp(-lambda))/factorial(x(i));
% lambda = mean / variance / E(x)
% lambda_optimal should = mean or sum(x)/n i.e. 6 in our data set
lambda = [0:12];
x = [12,4,5,3,7,5,6];
l_lambda_x=zeros(length(x),length(lambda));
for i =1:length(data_set)
l_lambda_x(i,:) =(x(i)*log(lambda))-log(factorial(x(i)))-lambda;
end
figure('visible','on')
plot(lambda,sum(l_lambda_x)');
[ll, pos]=max(sum(l_lambda_x));
lambda_hat = lambda(pos);


