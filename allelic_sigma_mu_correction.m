function [mu_minor, mu_major, sigma_mu, f]=allelic_sigma_mu_correction(hets,tau)
%error model for allelic capseg hets example of running function below


%model hets as sum of two beta distributions
%allele skew is corrected for in the tau already 
af_low=[0:.01:(.5)];
af_high=fliplr([.5:.01:1]);

%unused stats... 
hets.tumor_f=hets.i_t_alt_count./(hets.i_t_alt_count+hets.i_t_ref_count);
hets.shift=(abs(hets.tumor_f-.5));

p_down=1-poisscdf(sum(hets.tumor_f<.5),sum(hets.tumor_f>.5));
if sum(hets.tumor_f>.5)>0

while p_down<.3
    [v l]=min(hets.tumor_f);
    x=ones(slength(hets),1);
    x(l)=0;
    hets=reorder_struct(hets,x==1);
    p_down=1-poisscdf(sum(hets.tumor_f<.5),sum(hets.tumor_f>.5));
end

else
    mu_minor=NaN;
mu_major=NaN;
sigma_mu=NaN;
f=NaN;
return

end

for i=1:slength(hets)
    b_low(i,:)=betapdf(af_low,hets.i_t_alt_count(i)+1,hets.i_t_ref_count(i)+1);
    b_high(i,:)=betapdf(af_high,hets.i_t_alt_count(i)+1,hets.i_t_ref_count(i)+1);
end

%joint probability of each heterozygous SNP in the segment
b_log_sum=sum(log(b_high'+b_low'),2);
b_log_sum=exp(b_log_sum+(1-max(b_log_sum)));

b_log_sum=b_log_sum./sum(b_log_sum);
%lazy coding
beta_total=b_log_sum;

[~, peak_af]=max(beta_total);
sum_betas=b_high'+b_low';

outliers=repmat(.005,51,1);
for i=1:size(sum_betas,2)
    sum_betas(:,i)=sum_betas(:,i)./sum(sum_betas(:,i));
    if sum_betas(peak_af,i)>outliers(peak_af)
        out_hets(i)=0;
    else
        out_hets(i)=1;
    end
end
hets=reorder_struct(hets,out_hets==0);
    %try replacing outlier model with poission cdf # upper hets = lambda
    %(lambda fit?) weight the lowers by the distributions of lambdas
    %intergrate over all possible lambdas
    %#lower hets = x remove lower hets until p>.05~ or so. 

clear b_low b_high out_hets
if slength(hets)>1
for i=1:slength(hets)
    b_low(i,:)=betapdf(af_low,hets.i_t_alt_count(i)+1,hets.i_t_ref_count(i)+1);
    b_high(i,:)=betapdf(af_high,hets.i_t_alt_count(i)+1,hets.i_t_ref_count(i)+1);
end
b_log_sum=sum(log(b_high'+b_low'),2);
b_log_sum=exp(b_log_sum+(1-max(b_log_sum)));

b_log_sum=b_log_sum./sum(b_log_sum);
%lazy coding
beta_total=b_log_sum;

[~, peak_af]=max(beta_total);


int_beta_total=cumsum(beta_total);
m_p_up=min([1-int_beta_total(peak_af) .16]);
m_p_down=min([int_beta_total(peak_af) .16]);
%ensure that 1 sigma is captured regardless of peak location
ci_low=af_low(max([find(int_beta_total<(.32-m_p_up),1,'last') 1]));
ci_high=af_low(max([find(int_beta_total>.68+m_p_down,1,'first') peak_af]));


%overwriting allelic capseg fields
mu_minor=af_low(peak_af)*tau;
mu_major=(1-af_low(peak_af))*tau;
sigma_mu=((ci_high-ci_low))*tau;
f=af_low(peak_af);

else
    
mu_minor=NaN;
mu_major=NaN;
sigma_mu=NaN;
f=NaN;
end

end


function test
tumor_cov=load_table('~/Downloads/Tumor_0011.cov');
seg_file=load_table('~/Downloads/CLL-GCLL-0011-Tumor-SM-41JMA.tsv');
tumor_cov.Chromosome=chromosome2num_legacy(tumor_cov.Chromosome);

for i=1:slength(seg_file)
    hets=reorder_struct(tumor_cov,tumor_cov.Chromosome==seg_file.Chromosome(i)&tumor_cov.Start_position>=seg_file.Start_bp(i)&tumor_cov.Start_position<=seg_file.End_bp(i));
    if slength(hets)>1
    [seg_file.mu_minor(i,1), seg_file.mu_major(i,1), seg_file.sigma_major(i,1), seg_file.f(i,1)]=allelic_sigma_mu_correction(hets,seg_file.tau(i));
    seg_file.sigma_minor(i,1)=seg_file.sigma_major(i,1);
    else
        seg_file.mu_major(i,1)=NaN;
        seg_file.mu_minor(i,1)=NaN;
        seg_file.sigma_major(i,1)=NaN;
        seg_file.sigma_minor(i,1)=NaN;
    end
        
end
save_struct(seg_file,'CW103.seg')
end