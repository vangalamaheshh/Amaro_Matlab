
failed_genes = cell(size(fail_matrix, 1),1); 
for i = 1:size(fail_matrix,1) 
%    keyboard
    failed_genes{i} = M{1}.cov.targ.gene(find(M{1}.cov.targ.start == fail_matrix(i, 2) & M{1}.cov.targ.chr == fail_matrix(i, 1)));
end