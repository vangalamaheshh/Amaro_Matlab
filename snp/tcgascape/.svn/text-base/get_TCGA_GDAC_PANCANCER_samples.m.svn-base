function usamples = get_TCGA_GDAC_PANCANCER_samples
% get pan cancer samples from TCGA-GDAC 
segfile = ['/xchip/cga/gdac-prod/tcga-gdac/jobResults/GDAC_MergeDataFilesPipeline/PANCANCER/1035668/' ...
           '1.GDAC_MergeDataFiles.Finished/PANCANCER.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_cna__seg.seg.txt'];

fid = fopen(segfile);
usamples = [];
if fid > 0
    samples = textscan(fid,'%q %*[^\n]','TreatAsEmpty','NA',...
                       'Headerlines',1,'Delimiter',char(9)); 
    usamples = unique(samples{1}); 
end
fclose(fid);
