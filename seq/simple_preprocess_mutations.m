function M = simple_preprocess_mutations(M);

M.mut.patient = regexprep(M.mut.Tumor_Sample_Barcode,'-Tumor$','');
M.pat=[]; [M.pat.name tmp M.mut.pat_idx] = unique(M.mut.patient);
M.mut = add_simple_fieldnames(M.mut);
M.mut.newbase = find_newbase(M.mut);
M.mut.gene_idx = listmap(M.mut.gene,M.gene.name);
M.mut.context65 = get_from_fwb(M.mut.chr,M.mut.start,'/xchip/cga1/lawrence/db/hg19/context65/all.fwb');
[M.mut.categ M.mut.categ_ignoring_null_categ] = assign_mut_categs(M.mut,M.cat);

M.mut.is_flank = false(slength(M.mut),1); M.mut.is_flank(grepi('flank|utr|promoter|intron',M.mut.type,1)) = true;
M.mut.is_coding = false(slength(M.mut),1);
M.mut.is_coding(grepi('missense|nonsense|silent|splice|synon|frame|shift|non.?stop',M.mut.type,1)) = true;
M.mut.is_silent = false(slength(M.mut),1); M.mut.is_silent(grepi('silent|synon',M.mut.type,1)) = true;

