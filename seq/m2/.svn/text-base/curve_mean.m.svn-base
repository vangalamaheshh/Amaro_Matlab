function [ mean_curve ] = curve_mean(length)

%    mean_curve = zeros(length, 1);
%    for i = 1:length
%        mean_curve(i) = npairs*(2*(i/length) - (i/length)^2); 
%    end
 
a = (1:length)'/length;
mean_curve = 2*a - a.^2;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% code to test it

n=10000;
x=ceil(length*rand(n,1));
d=abs(bsxfun(@minus,x,x'));
h=histc(d(:),0:(length-1));
h(1)=h(1)-n; % subtract diagonal
c=cumsum(h)/sum(h);
pr(mean_curve,c,1:length)


% RESULT: approaches a perfect match as length increases
