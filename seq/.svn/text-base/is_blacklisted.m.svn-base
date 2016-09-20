function isb = is_blacklisted(maf,varargin)
% is_blacklisted(maf,varargin)
% use the same as:
% apply_mutation_blacklist(maf,varargin)

if ischar(maf), maf=load_struct(maf); end

maf.temporary_orig_maf_idx = (1:slength(maf))';
maf2 = apply_mutation_blacklist(maf,varargin{:});

isb = ~ismember((1:slength(maf))',maf2.temporary_orig_maf_idx);

