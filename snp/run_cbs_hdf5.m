function run_cbs_hdf5(prefix,hdf5file,chrposfile,samples_idx,paramsstr,sampleID)

if ~exist('sampleID','var') || isempty(sampleID)
    sampleID = prefix;
end

fid = fopen([prefix '.R'],'w');
fprintf(fid,'source("~/CancerGenomeAnalysis/trunk/R/getCBSparam_hdf5.R")\n');
fprintf(fid,['getCBSparam_hdf5("' chrposfile '","' prefix '.seg.dat", "' ...
    hdf5file '","/dat",' num2str(samples_idx) ',"' sampleID '",' paramsstr ')' newline]);
fprintf(fid,'q()');

fclose(fid);