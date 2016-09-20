function n=sample_size(fhandle,p,alpha,p0,tail,power_thresh)
%x0=fzero(@(x) norminv(power_thresh,exp(x)*p,...
%                      sqrt(exp(x)*p*(1-p)))-...
%         norminv(alpha,exp(x)*p0,...
%                 sqrt(exp(x)*p0*(1-p0))),log(50));
x=fzero(@(x) find_N(fhandle,p,alpha,x,p0,tail,power_thresh),log(50)); %was x0 
n=floor(exp(x));

function res=find_N(fhandle,p,alpha,lgN,p0,tail,power_thresh)
[res,res1]=feval(fhandle,p,alpha,p0,floor(exp(lgN)),tail);
res=res-power_thresh;

function res=find_approx_N(fhandle,p,alpha,lgN,p0,tail,power_thresh)
res=feval(fapprox_handle,p,alpha,p0,floor(exp(lgN)),tail)-power_thresh;
