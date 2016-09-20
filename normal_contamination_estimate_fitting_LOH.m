function [log_odds, p2, c1]=normal_contamination_estimate_fitting_LOH(t_alt_count,t_ref_count,n_alt_count,n_ref_count,fit1,fit2,dx)
    af_range=0:dx:1; 
     
    normalf=n_alt_count/(n_alt_count+n_ref_count);
    tumorf=t_alt_count/(t_alt_count+t_ref_count);
    
    
   model1=af_range;
   
%    if (normalf<.5 && tumorf>.5) || (normalf>.5 && tumorf<.5)
%     model2=(af_range-.5)*fit2+normalf;
%    else
       model2=(af_range-.5)*fit2+.5;
%    end
   
    %downsample so that the tumor and the normal have the same coverage
    normcov=n_alt_count+n_ref_count;
    tumcov=t_alt_count+t_ref_count;
    
    % gaussian noise around .48 std .02 .02*randn(100,1)+.48
    
    
    % Probability of a give allele fraction given t_alt_count
    % add 2 parameters (0,.02) fitting beta distribution mean 0 dist of .02
    wt=betapdf(af_range,t_alt_count+1,t_ref_count+1)*dx;
    %wnt=betapdf(af_range,
    % Squished by TiN to project on to normal allele fractions
    %wt=conv(wt,normpdf(0:.001:1,.48,.04),'same');
    
    hwt=histw(wt,model2,model1); %weights for the normal projection
   %     bwn=binopdf(n_alt_count,normcov,af_range);
   % bwn=bwn/sum(bwn);
   % p4=sum(hwt.*bwn);
    
    %wn2=betapdf(model2,n_alt_count+1,n_ref_count+1);
    
    % Probability of a given allele fraction give n_alt_count 
    
    
    wn1=betapdf(model1,n_alt_count+1,n_ref_count+1)*dx;
   % wn1=conv(wn1,normpdf(0:.001:1,.48,.04),'same');
  
    
   
    c1=sum(wn1);
    c2=sum(hwt);
    c3=sum(wt);
    if c1<.95||c2<.95||c3<.95
    disp('Probabilities not summing up to 1!')
    end
    
    p1=wn1(round(length(af_range)/2));
    p2=dot(hwt,wn1);
    
    
    
    log_odds=log(p1)-log(p2);
    
end
function test


normal_contamination_estimate_fitting_LOH(46,90,76,79,1,0,.001)

end


