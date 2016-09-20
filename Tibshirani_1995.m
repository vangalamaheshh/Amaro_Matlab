%Regression Shrinkage and Selection via Lasso
%y response x i1-ip predictors

load('longely.mat');
x=longely(:,1:6);
y=longely(:,7);
N=16;
% Beta1...p coefficents
% lamda tuning paramter
% standardize x
% alpha constant

for j=1:size(x,2)
    x(:,j)=(x(:,j)-mean(x(:,j)))/std(x(:,j));
end

% (a,B) = arg min( for all i...N sum(yi - a - sum over columns Bjxij)^2 )
% subject to |Bj|<=t
% a = mean(y) mean(y)=0 therefore omit a

%1/2N * 
alpha=0;
for i=1:N
 least_sq=(y(i)-x(i,:)*Beta)^2;
end
pen=least_sq+lambda*sum(abs(Beta));
