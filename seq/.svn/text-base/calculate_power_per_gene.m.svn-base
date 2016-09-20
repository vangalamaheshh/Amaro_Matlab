function calculate_power_per_gene(M,P)
% based on code by Gaddy Getz, 2010

require_fields(M,{'gene','ng','patient','np','N_terr','TOT','mutrate'});

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'power_per_gene_outstem','*required*');
P = impose_default_value(P,'patient_fracs_to_test',[0.03 0.05 0.1 0.15 0.2]);
P = impose_default_value(P,'FNR',0.15);

% addpath ~gadgetz/projects/power     % for single_screen_power_smooth

fprintf('Analyzing power per gene...\n');

tot_outfile = [P.power_per_gene_outstem '.total_coverage.txt'];
hist_outfile = [P.power_per_gene_outstem '.total_coverage_histogram'];
pow_outfile = [P.power_per_gene_outstem '.power_stats.txt'];
powfull_outfile = [P.power_per_gene_outstem '.power_stats_full.txt'];
summary_outfile = [P.power_per_gene_outstem '.power_summary.txt'];

if isfield(M,'N_cov')
  actual_cov = squeeze(M.N_cov(:,M.TOT,:));
else  
  fprintf('Imputing full coverage in calculate_power_per_gene\n');
  actual_cov = squeeze(M.N_terr(:,M.TOT,:));
end
total_cov = sum(actual_cov,1);
f=fopen(tot_outfile,'w');
for i=1:length(total_cov), fprintf(f,'%s\t%ld\n',M.patient.name{i},round(total_cov(i))); end
fclose(f);

figure(1); clf;
hist(total_cov/1e6,50);
title('total sequencing coverage');
xlabel('Coverage [Mbp]');
ylabel('Count');
print('-dpng','-r180',hist_outfile);
%print_D(hist_outfile,{{'pdf'},{'png','-r180'}});  % (was causing problems in Firehose)

frac_of_patients=P.patient_fracs_to_test;
false_neg_rate=P.FNR;

BMR = M.mutrate.tot.hat;
mu=BMR;
alpha=0.1/M.ng;

for full=0:1
 if full, f=fopen(powfull_outfile,'w');
 else f=fopen(pow_outfile,'w');
 end
 fprintf(f,'Gene\tLength\t# of samples\tTotal territory\tCovered territory\tFalse negative rate\tMissing data rate\tBMR\tAlpha');
 for j=1:length(frac_of_patients), fprintf(f,'\t%d',round(frac_of_patients(j)*100)); end, fprintf(f,'\n');
  for i=1:M.ng
    len = M.N_terr(i,M.TOT);
    total = len * M.np;   
    covered = sum(actual_cov(i,:));
    if full, missing_rate=false_neg_rate;
    else missing_rate=1-(covered/total)*(1-false_neg_rate);
    end
    p=zeros(1,length(frac_of_patients));
    for k=1:length(frac_of_patients)
      try
        p(k)=single_screen_power_smooth(mu*len,mu*len+(1-missing_rate)*frac_of_patients(k),M.np,M.ng,1,0.1,len,150);
      catch me
        fprintf('ERROR in single_screen_power_smooth for gene %s:\n%s\n',M.gene.name{i},me.message);
        p(k) = nan;
      end
    end
    fprintf(f,'%s\t%ld\t%ld\t%ld\t%ld\t%f\t%f\t%f\t%f',M.gene.name{i},len,M.np,total,covered,false_neg_rate,missing_rate,BMR,alpha);
    for j=1:length(p), fprintf(f,'\t%f',p(j));end, fprintf(f,'\n');
    if mod(i,1000)==0, fprintf('%d/%d ',i,M.ng); end
  end,fprintf('\n');
 fclose(f);
end

% summary of how many genes we have 80% power for

x1 = load_struct(pow_outfile);
x2 = load_struct(powfull_outfile);
x1 = make_numeric(x1,{'x3','x5','x10','x15','x20'});
x2 = make_numeric(x2,{'x3','x5','x10','x15','x20'});
p1 = [x1.x3 x1.x5 x1.x10 x1.x15 x1.x20];
p2 = [x2.x3 x2.x5 x2.x10 x2.x15 x2.x20];
pow = [3 5 10 15 20];
f = fopen(summary_outfile,'wt');
fprintf(f,'%% patients w/mut\tActual coverage\t\tFull coverage\n');
for i=1:length(pow)
  fprintf(f,'%d%%\t\t%d/%d\t\t%d/%d\n',pow(i),...
          sum(p1(:,i)>0.8),size(p1,1),...
          sum(p2(:,i)>0.8),size(p2,1));
end
fclose(f);


