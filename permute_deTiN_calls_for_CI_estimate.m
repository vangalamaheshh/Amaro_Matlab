function [  sx ] = permute_deTiN_calls_for_CI_estimate(call,k )
%calculate the point estimate while excluding random mutations used in the
%fit; store these fits and construct CI's based on permuted point estimates
perm=reorder_struct(call,k==1);
disp('performing bootstrap confidence interval estimation')
rng(1)
perm_TiN=zeros(100,1);
%%
% Sim_for_CRSP=[.01,.05,.1,.2,.3,.5,.8];
% for j=1:length(Sim_for_CRSP)
for i=1:100
    index=zeros(slength(perm),1);
    index(randi(slength(perm),round(slength(perm)*.8),1))=1;
    [perm_TiN(i)]=loglikelihood_fit(perm,index==1);
end
% j
% end
TiN=mean(perm_TiN);
sx=std(perm_TiN);
if TiN-(2*sx)>0.01
    CI(1)=TiN-(2*sx);
else
    CI(1)=0;
end
if TiN+(2*sx)>1
    CI(2)=1;
else
    
    CI(2)=TiN+(sx*2);
    
end



end

