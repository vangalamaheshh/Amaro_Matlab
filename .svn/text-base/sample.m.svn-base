function s=sample(n,k)
% sample k without replacement from n
% if n is a number then assumes the set is 1:n

if length(n)==1
  n=1:n;
end

r=randperm(length(n));
s=n(r(1:k));


