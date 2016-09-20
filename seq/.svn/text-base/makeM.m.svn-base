function M = makeM(m)

M=[];
M.mut = add_and_convert_simple_fieldnames(m);
M.mut = add_helper_is_fields(M.mut);

[M.pat.name ui M.mut.pat_idx] = unique(M.mut.patient);
M.pat.nmut = histc(M.mut.pat_idx,1:slength(M.pat));

if isfield(M.mut,'ttype')
  M.pat.ttype = M.mut.ttype(ui);
  [M.ttype.name vi M.pat.ttype_idx] = unique(M.pat.ttype);
  M.mut.ttype_idx = M.pat.ttype_idx(M.mut.pat_idx); 
  M.ttype.nmut = histc(M.mut.ttype_idx,1:slength(M.ttype));
end

