
lift_over_hg19 = fopen('lift_over_hg19.txt', 'wt');
for i = 1:size(M{1}.cov.targ.gene,1)
    if M{1}.cov.targ.chr(i) ~= 24 & M{1}.cov.targ.chr(i) ~= 23
        fprintf(lift_over_hg19, 'chr%d:%d-%d\n', M{1}.cov.targ.chr(i), M{1}.cov.targ.start(i), M{1}.cov.targ.end(i));
    elseif M{1}.cov.targ.chr(i) == 24
        fprintf(lift_over_hg19, 'chrY:%d-%d\n', M{1}.cov.targ.start(i), M{1}.cov.targ.end(i));
    else 
        fprintf(lift_over_hg19, 'chrX:%d-%d\n',  M{1}.cov.targ.start(i), M{1}.cov.targ.end(i));
    end
end