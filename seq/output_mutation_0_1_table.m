function X = output_mutation_0_1_table(M,filename)

X=[];
X.gene = M.gene.name;
if isfield(M.gene,'chr'), X.chr = M.gene.chr; end
if isfield(M.gene,'start'), X.start = M.gene.start; end
if isfield(M.gene,'end'), X.end = M.gene.end; end
if isfield(M.gene,'len'), X.len = M.gene.len; end

pname = regexprep(M.patient.name,'-','_');
pname = genfieldname(pname);

for i=1:M.np
  X = setfield(X,pname{i},double(M.n_nonsilent(:,M.TOT,i)>0));
end

save_struct(X,filename);

