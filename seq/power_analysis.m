function p=power_analysis(mu,len,n_samples,n_genes,missing_data_rate,frac_of_patients_w_mut)

if (ischar(mu) || isstruct(mu)), test_power_analysis(mu); return; end

p=single_screen_power_smooth(mu*len,mu*len+(1-missing_data_rate)*frac_of_patients_w_mut,...
                             n_samples,n_genes,1,0.1,len,150);


%%%%%%%
%% test power_analysis

function test_power_analysis(params)
if ischar(params)
  tmp=struct('method',params);
  params=tmp;
end

switch params.method
 case 'simple'
  p=power_analysis(1.5e-6,1500,500,20000,0.23,0.03); % 88.97%
  
 case 'mcfadden'
  for mu=[0.5 1]*1e-6
    for ng=[2000 20000]
      for ns=[25 50]
        for sig=[0.05 0.1 0.15 0.2]
          p=power_analysis(mu,1500,ns,ng,0.15,sig); % 
          fprintf(1,'mu=%5.3f [mut/Mb]\tn_samples=%d\tn_genes=%d\tsignal=%f\t==>\tpower=%f\n',mu*1e6,ns,ng,sig,p);
        end
      end
    end
  end
  
 otherwise 
  error('no such method');
end
