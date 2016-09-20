function M = maftoM(maf)

if ischar(maf), maf = load_struct(maf); end

demand_fields(maf,{'patient','ttype'});

M=[]; M.mut=maf;
[M.pat.name ui M.mut.pat_idx] = unique(M.mut.patient);
M.pat.ttype = M.mut.ttype(ui);
M.pat.nmut = histc(M.mut.pat_idx,1:slength(M.pat));
[M.ttype.name ui M.pat.ttype_idx] = unique(M.pat.ttype);
M.ttype.npat = histc(M.pat.ttype_idx,1:slength(M.ttype));
for i=1:slength(M.ttype), M.ttype.nmut(i,1) = sum(M.pat.nmut(M.pat.ttype_idx==i)); end

