function c = load_coverage1040(sample,params)
if iscell(sample), error('Multiple samples not supported'); end

if ~exist('params','var'), params=[]; end
params = impose_default_value(params,'sum_chromosomes',true);

fn = ['/xchip/tcga_scratch/lawrence/' sample '/genome_allcateg_coverage.txt'];
demand_file(fn);
X = read_table(fn,['%s' repmat('%f',1,3+1040)],char(9),0,'whitespace','\b\r');

c = cat(2,X.dat{5:end});

if params.sum_chromosomes
  c = sum(c,1);
end

c = as_column(c);
