function [CI_l, CI_h,curve] =  permute_deTiN_combined_estimate(call,k,b,call_filtered,cc_TiN)
permutations = 100;
b=rmfield_if_exist(b,{'header','headline'});
call = reorder_struct(call,k);
TiN=0:.01:1;
b = reorder_struct(b,b.clust_assignment==mode(b.clust_assignment));
b.expected_shift=zeros(slength(b),2,length(TiN));
call_filtered.direction = str2double(call_filtered.direction);
call_filtered.x=str2double(call_filtered.x);
call_filtered.seg_id=zeros(slength(call_filtered),1);



for seg=1:slength(b)
    cfilt=reorder_struct(call_filtered,call_filtered.x>=b.xs(seg)&call_filtered.x<=b.xe(seg));
    cfilt.direction=cfilt.tumor_f>mean(cfilt.normal_f);
    AFt(1)=mean(cfilt.tumor_f(cfilt.direction==1))-mean(cfilt.tumor_f);
    AFt(2)=mean(cfilt.tumor_f)-mean(cfilt.tumor_f(cfilt.direction==0));
    CNt = b.tau(seg); CNn = 2;
    b.expected_shift(seg,1,:)=AFt(1)*(CNt*TiN)./(((1-TiN)*CNn)+(TiN*CNt));
    b.expected_shift(seg,2,:)=AFt(2)*(CNt*TiN)./(((1-TiN)*CNn)+(TiN*CNt));
    call_filtered.seg_id(call_filtered.x>=b.xs(seg)&call_filtered.x<=b.xe(seg),1)=seg;
end
call_filtered = reorder_struct(call_filtered,~call_filtered.seg_id==0);


events = slength(call)+slength(call_filtered);
rng(1)
perm_pTiN=zeros(permutations,length(TiN));
mutation_likelihood=zeros(permutations,length(TiN));
for p_i=1:permutations
    index_SSNVs=zeros(slength(call),1);
    index_SSNVs(randi(slength(call),round(slength(call)*.95),1))=1;
    index_SNPs=zeros(slength(call_filtered),1);
    index_SNPs(randi(slength(call_filtered),round(slength(call_filtered)*.95),1))=1;
    
    perm=reorder_struct(call_filtered,index_SNPs==1);
    pTiN=zeros(slength(perm),length(TiN));
    for i=1:slength(perm)
        if perm.direction(i)==1
            pTiN(i,:)=betapdf(mean(perm.normal_f)+b.expected_shift(perm.seg_id(i),1,:),perm.n_alt_count(i)+1,perm.n_ref_count(i)+1);
        else
            pTiN(i,:)=betapdf(mean(perm.normal_f)-b.expected_shift(perm.seg_id(i),2,:),perm.n_alt_count(i)+1,perm.n_ref_count(i)+1);
        end
        
        pTiN(i,:)=pTiN(i,:)./sum(pTiN(i,:));
        
        
    end
    pTiN=sum(log(pTiN));
    pTiN=pTiN+(1-max(pTiN));
    pTiN=exp(pTiN);
    perm_pTiN(p_i,:)=pTiN;
    
    
    [perm_TiN(i),mutation_likelihood(p_i,:)]=loglikelihood_fit(call,index_SSNVs==1);
    
    
end
combined_likelihood=log(perm_pTiN')-log(mutation_likelihood');
[v arg] = max(combined_likelihood);
TiN_est =  cc_TiN;
CI_l = max([0, TiN_est-(2*std(TiN(arg)))]);
CI_h = min([1, TiN_est+(2*std(TiN(arg)))]);
curve = sum(combined_likelihood,2);
end




