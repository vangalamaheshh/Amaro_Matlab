function [wa wb] =mle_beta_w(x,X,n,N)
%  (a,b) --> # r --> [[n]] <-- [[N]] #_k
% log[ L(a,b|{n,N}_1..k) ] = 
% Sum_k{ log [ Int( Binom(n|N,r) Beta(r|a,b) dr ) ] } =
% Sum_k{ log [ hype2pdf(n,N,a,b) ] } = 
% Sum_k{ hyge2ln(n,N,a,b) }

w1=0:0.01:1;
l=nan(length(w1),1);
for i=1:length(l)
  l(i)=bbinologlik(x,X,w1(i)*n+1,w1(i)*(N-n)+1);
end

plot(w1,l);

[tmp idx] = max(l);
wa = w1(idx);

dist_params=0.3;

fprintf('Initial guess:  w=%f  bbinologlk=%f\n',...
        dist_params(1), -bbinologlik(x,X,dist_params*n+1,dist_params*(N-n)+1));

xmle=fminsearch( @(x) -bbinologlik(x,X,dist_params*n+1,dist_params*(N-n)+1),dist_params,optimset('Display','iter','MaxFunEvals',2000));

wb=xmle;



