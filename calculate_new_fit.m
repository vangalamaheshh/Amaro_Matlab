function [fit_n,k,gof] = calculate_new_fit (a,call,germline_fit,threshold,kold)

germline_fit=1;

% log ratio consistent with point?
for f=1:slength(call)
call.iter(f,1)=normal_contamination_estimate_fitting(call.t_alt_count(f),call.t_ref_count(f)...
         ,call.n_alt_count(f),call.n_ref_count(f),a,germline_fit);

end

sites_to_consider=(ismember(call.judgement,'KEEP')|...
         ismember(call.failure_reasons,'germline_risk,normal_lod,alt_allele_in_normal')| ...
          ismember(call.failure_reasons,'normal_lod,alt_allele_in_normal') | ismember(call.failure_reasons,'normal_lod') | ismember(call.failure_reasons,'alt_allele_in_normal') |...
          ismember(call.failure_reasons,'germline_risk') | ismember(call.failure_reasons,'alt_allele_in_normal,strand_artifact'));
     
k=(call.iter>threshold)&sites_to_consider;

g = fittype('a*x');
[fit_n,gof]=fit(call.tumor_f(k),call.normal_f(k),g,'Weights',call.w_beta(k));





end


function test
sg=1;
fit_n.a=outfit.a;
i=1;
slopes(1,1)=outfit.a;
N_k(1,1)=sum(k);
max_iter=100;
f_c=10^-4;
while sg
    i=i+1;
    [fit_n,k]=calculate_new_fit(fit_n.a,call,linfit.a,1.5,k);
    slopes(i,1)=fit_n.a;
    N_k(i,1)=sum(k);
    delta=(slopes(i)-slopes(i-1))/mean(slopes((i-1):i));
    sg=(i < max_iter)&(f_c<delta);
    slopes(i)

end

end