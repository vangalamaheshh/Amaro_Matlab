function [a,b]=mle_beta(n,N)
%  (a,b) --> # r --> [[n]] <-- [[N]] #_k
% log[ L(a,b|{n,N}_1..k) ] = 
% Sum_k{ log [ Int( Binom(n|N,r) Beta(r|a,b) dr ) ] } =
% Sum_k{ log [ hype2pdf(n,N,a,b) ] } = 
% Sum_k{ hyge2ln(n,N,a,b) }

if (0)
    r=n./N;
    r(isnan(r) | isinf(r))=[];
    dist_params=betafit(r);
    dist_params
end

dist_params=[2 2*(sum(N)/sum(n)-1)];

fprintf('Initial guess:  a=%f  b=%f  bbinologlk=%f\n',...
        dist_params(1),dist_params(2), -bbinologlik(n,N,dist_params(1),dist_params(2)));

% xmle=fminsearch( @(x) -hyge2loglik(n,N,x(1),x(2)),dist_params,optimset('Display','iter','MaxFunEvals',2000));
xmle=fminsearch( @(x) -bbinologlik(n,N,x(1),x(2)),dist_params,optimset('Display','iter','MaxFunEvals',1000));

%xmle=fminsearch( @(x) -hyge2loglik(n,N,x(1),x(2)),dist_params,optimset('MaxFunEvals',2000));

a=xmle(1);
b=xmle(2);


