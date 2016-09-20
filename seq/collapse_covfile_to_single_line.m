function collapse_covfile_to_single_line(samples,covfile)
% collapse_covfile_to_single_line(samples,covfile)
%
% Mike Lawrence 2010-01-27

cmds={}; banners={};
for i=1:length(samples)
  in = ['/xchip/tcga_scratch/lawrence/' samples{i} '/' covfile];
  out = regexprep(in,'\.txt$','\.tot.txt');
  if strcmp(in,out), error('"covfile" should end in ".txt"'); end
  if ~exist(out,'file')
    cmds{end+1} = ['"more ' in ' | totcols > ' out '"'];
    banners{end+1} = ['TOTCOLS' num2str(i)];
  end
end
bwait(bsub(cmds,banners));
