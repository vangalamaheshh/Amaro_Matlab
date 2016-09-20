function [pmax pmid]=truncated_product_combined_p(pvec,tau)

if ~exist('tau','var'), tau=0.05; end

if numel(pvec)==1, error('needs at least two p-values'); end
if size(pvec,2)==1 && size(pvec,1)>1, pvec=pvec'; end

if size(pvec,1)>1
  pmax = nan(size(pvec,1),1);
  pmid = nan(size(pvec,1),1);
  for i=1:size(pvec,1)
    [pmax(i) pmid(i)] = truncated_product_combined_p(pvec(i,:),tau);
  end
  return
end

n = length(pvec);
pvec(pvec>=tau)=1;
sc=-log(pvec);
sc_obs=sum(sc);

if sc_obs==0
  pmax = 1; 
  q = (1-tau)^2;
  pmid = q + (1-q)*rand;
  return
end

method=2;

if method==1

  % MANUAL CONVOLUTIONS METHOD

  step_size=1e-8;
  pnull=step_size:step_size:tau;

  pnull(pnull>=tau)=1;
  sc_null=-log(pnull);
  
  sub_bins=1000;
  sc_step=abs(min(diff(sc_null(1:(end-1)))))/sub_bins;
  sc_hist=histc(sc_null,0:sc_step:max(sc_null))/(length(pnull)/tau);
  sc_hist(1)=1-tau;
  
  sc_conv=sc_hist;
  for i=2:n
    sc_conv=conv(sc_conv,sc_hist);
  end
  %disp(length(sc_conv))
  %disp(sc_step)
  %disp(sc_obs)

  pidx=floor(sc_obs/sc_step);  % make sure that the indices are not off by 1
  pmax=1-sum(sc_conv(1:pidx));
  
elseif method==2

  % ANALYTICAL METHOD

  % for length of two
  % single func:
  % 1-tau | x=0   
  % exp(-x) | a<=x<Inf

  % conv with itself:
  % (1-tau)^2 |  x=0
  % 2*(1-tau)*exp(-x) | a<=x< 2*a
  % (x-2*a)*exp(-x)+ 2*(1-tau)*exp(-x)   | 2*a <= x < Inf

  if n~=2 && n~=3, error('not supported'); end

  pmax=1-tpmcdf(sc_obs,n,tau);

else
  error('invalid method');
end

pmid = pmax;
