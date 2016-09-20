function [sigma] =  permute_deTiN_LOH_snps_seg_estimate(cfilt,expected_shift,seg_id)

rng('default')
TiN=0:.01:1;
cfilt=reorder_struct(cfilt,ismember(cfilt.seg_id,seg_id));
perm_pTiN=zeros(100,101);

for p_i=1:100
    
    index=zeros(slength(cfilt),1);
    index(randi(slength(cfilt),round(slength(cfilt)*.95),1))=1;
    perm=reorder_struct(cfilt,index==1);
    
    pTiN=zeros(slength(perm),101);
    for i=1:slength(perm)
        if perm.direction(i)==1
            pTiN(i,:)=betapdf(mean(perm.normal_f)+expected_shift(seg_id==perm.seg_id(i),1,:),perm.n_alt_count(i)+1,perm.n_ref_count(i)+1);
        else
            pTiN(i,:)=betapdf(mean(perm.normal_f)-expected_shift(seg_id==perm.seg_id(i),2,:),perm.n_alt_count(i)+1,perm.n_ref_count(i)+1);
        end
        
        pTiN(i,:)=pTiN(i,:)./sum(pTiN(i,:));
        
        
    end
    pTiN=sum(log(pTiN));
    pTiN=pTiN+(1-max(pTiN));
    pTiN=exp(pTiN);
    perm_pTiN(p_i,:)=pTiN;
    
end

[value perm_TiN] = max(perm_pTiN');
perm_TiN=TiN(perm_TiN);
CI_low = max([0,(mean(perm_TiN)-(2*std(perm_TiN)))]);
CI_high = min([1,(mean(perm_TiN)+(2*std(perm_TiN)))]);
sigma = std(perm_TiN);
TiN_est=mean(perm_TiN);
end