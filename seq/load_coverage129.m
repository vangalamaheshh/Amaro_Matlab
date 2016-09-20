function c = load_coverage129(sample)
if iscell(sample), error('Multiple samples not supported'); end

fn = ['/xchip/tcga_scratch/lawrence/' sample '/genome_context129_coverage.txt'];
demand_file(fn);
X = read_table(fn,['%s' repmat('%f',1,3+129)],char(9),0,'whitespace','\b\r');

c = cat(2,X.dat{5:end});

