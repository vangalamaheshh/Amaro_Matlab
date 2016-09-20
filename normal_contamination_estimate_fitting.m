function [log_odds, pTiN, pGerm]=normal_contamination_estimate_fitting(t_alt_count,t_ref_count,n_alt_count,n_ref_count,TiNfit,Germfit,dx,tau)
    
    af_range=(dx/2)-10e-20:dx:1-(dx/2)+10e-20; 
    tau=tau+.0001; %handle tau =0
    CNR=tau/2;
    TiNo=CNR/(CNR+(1/TiNfit)-1);
    afTiN=TiNo*af_range;
    afGerm=Germfit*af_range;
    n_depth=n_alt_count+n_ref_count;
    Ex_wt = 1/(1+(t_ref_count/t_alt_count));
    TiN_allele_count = Ex_wt*TiNfit*n_depth;
    %This is the allele fraction probabilities in the af range for both T
    %and N
    wt=betapdf(af_range,t_alt_count+1,t_ref_count+1)*dx;
    wn=betapdf(af_range,n_alt_count+1,n_ref_count+1)*dx;
    wTiN=betapdf(af_range,TiN_allele_count+1,n_depth-TiN_allele_count+1)*dx;
   % chopped proability solution by fit1 original -- Legacy
   % wn1=betapdf(model1,n_alt_count+1,n_ref_count+1);
    
    
    
   % correction for low allele count bias 
%     wx=binocdf(3,t_alt_count+t_ref_count,af_range);
%     wt=(wt.*(wx+1))./(sum(wt.*(wx+1)));
%     
    %projecting the tumor into TiN normal space
   % wTiN=betapdf(afGerm,
   % wTiN=histw(wt,afTiN,afGerm);
    

    wn=wn/sum(wn);
    wTiN=wTiN/sum(wTiN);
    wt=wt/sum(wt);
    
    
    pTiN=dot(wTiN,wn);
    

    pGerm=dot(wt,wn);
    
    
    
    log_odds=log(pTiN)-log(pGerm);
    
  
    
   
  
end
function test


normal_contamination_estimate_fitting(180,233,0,200,0,1,.01,2)

end


