function fh_RearrangementsMultiplicity(sample,tumor_cr_fn,normal_cr_fn,purity,ploidy,dRangerFile)
% Yotam Drier, yotamd@gmail.com

if nargin~=6, error('Usage: fh_RearrangementsMultiplicity(sample,tumor_cr_fn,normal_cr_fn,purity,ploidy,dRangerFile)'); end

demand_file(tumor_cr_fn);
demand_file(normal_cr_fn); 
if ~isnumeric(purity), purity=str2double(purity); end
if ~isnumeric(ploidy), ploidy=str2double(ploidy); end

P=[];
P.purity=purity;
P.ploidy=ploidy;
RearrangementsMultiplicity(sample,tumor_cr_fn,normal_cr_fn,P,dRangerFile) 


