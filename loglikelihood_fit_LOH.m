function [l, noise_subtracted_ll_curve, call]=loglikelihood_fit_LOH(call)
a=-1:.01:1;


a_log=zeros(length(a),1);


for j=1:length(a)
    if mod(a(j),.2)==0
        disp(sprintf('Checking TiN=%d',a(j)));
    end
    for i=1:slength(call)
    
    [ll_het_consistent_with_TiN p_TiN]=normal_contamination_estimate_fitting_LOH(call.t_alt_count(i),call.t_ref_count(i)...
            ,call.n_alt_count(i),call.n_ref_count(i),1,a(j),.001);
    a_log(j,1)=a_log(j,1)+ll_het_consistent_with_TiN;
    call.pTiNs(i,j)=p_TiN;
    call.ll_TiN(i,1)=ll_het_consistent_with_TiN;
    end
end
noise_subtracted_ll_curve=sum(call.pTiNs(:,101:end))-(sum(call.pTiNs(:,1:101)));
[v l]=max(noise_subtracted_ll_curve);
l=l-1;


end