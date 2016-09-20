function ci=binoci(x,n,alpha,initial_guess)
% Gaddy Getz 2010

if ~exist('alpha','var') || isempty(alpha)
    alpha=0.05;
end

if ~exist('initial_guess','var') || isempty(initial_guess)
    [tmp,initial_guess]=binofit(x,n,alpha);
    df=binoci_diff(initial_guess,x,n,alpha);
    disp(['binofit: ' num2str(initial_guess) ' ; interval=' num2str(abs(diff(initial_guess))) ' ; diff=' num2str(df)]);
end

ci(1)=norminv(alpha/2,x/n,sqrt((x/n)*(1-(x/n))/n));
ci(2)=norminv(1-alpha/2,x/n,sqrt((x/n)*(1-(x/n))/n));
df=binoci_diff(ci,x,n,alpha);
disp(['norminv: ' num2str(ci) ' ; interval=' num2str(abs(diff(ci))) ' ; diff=' num2str(df)]);

ci(1)=betainv(alpha/2,x+1,n-x+1);
ci(2)=betainv(1-alpha/2,x+1,n-x+1);
df=binoci_diff(ci,x,n,alpha);
disp(['betainv: ' num2str(ci) ' ; interval=' num2str(abs(diff(ci))) ' ; diff=' num2str(df)]);

ci=fsolve(@(t) binoci_diff(t,x,n,alpha),initial_guess,optimset('Display','off','TolFun',1e-9));
df=binoci_diff(ci,x,n,alpha);
disp(['min_beta_interval: ' num2str(ci) ' ; interval=' num2str(abs(diff(ci))) ' ; diff=' num2str(df)]);
