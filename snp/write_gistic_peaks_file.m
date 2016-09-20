function write_gistic_peaks_file(fname,peaks)

npeaks = length(peaks);
longest_peak = max(cellfun(@length,{peaks.genes}));
gene_grid = repmat({''},longest_peak,npeaks);
for p = 1:npeaks
    genes = peaks(p).genes;
    ngenes = length(genes);
    gene_grid(1:ngenes,p) = genes';
end
fid = fopen(fname,'w');
fprintf(fid,'cytoband');
fprintf(fid,'\t%s',peaks.cytoband);
fprintf(fid,'\n');

fprintf(fid,'q value');
fprintf(fid,'\t%g',peaks.qv);
fprintf(fid,'\n');

fprintf(fid,'residual q value');
fprintf(fid,'\t%g',peaks.resid_qv);
fprintf(fid,'\n');

fprintf(fid,'wide peak boundaries');
fprintf(fid,'\t%s',peaks.loctext);
fprintf(fid,'\n');

fprintf(fid,'genes in wide peak');
for i = 1:longest_peak
    fprintf(fid,'\t%s',gene_grid{i,:});
    fprintf(fid,'\n');
end
fclose(fid);