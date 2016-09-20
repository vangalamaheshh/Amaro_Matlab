function [t n] = get_bam_stats(sample)

try

basedir = ['/xchip/tcga_scratch/lawrence/' sample];
L = {'in total','QC failure','duplicates','mapped','paired','read1','read2','properly paired',...
  'with itself and mate mapped','singletons','with mate mapped to a different chr',...
  'with mate mapped to a different chr (mapQ>=5)'};
F = {'total','non_pf','duplicates','mapped','paired','read1','read2','properly_paired',...
  'both_mapped','singletons','interchr','interchr_mapQ5'};

if length(L) ~= length(F), error('get_bam_stats: L and F must be same length'); end

for i=1:2
  if i==1, tn = 'tumor'; else tn = 'normal'; end
  file = [basedir '/' tn '.stats'];
  if ~exist(file,'file')
    flagstat(sample);
  end
  x = load_lines(file);
  s = [];
  for j=1:length(x)
    spc = find(x{j}==' ',1);
    if isempty(spc), continue; end
    if spc==length(x{j}), continue; end
    n = str2double(x{j}(1:spc-1));
    m = x{j}(spc+1:end);
    match = [];
    for k=1:length(L)
      if strncmpi(m,L{k},length(L{k})), match=[match;k]; end
    end
    if isempty(match), continue; end
    if length(match)>1
      [tmp ord] = sort(cellfun('length',L(match)));
      match = match(ord(end));
    end
    % fprintf('%s -> %s\n',L{match},F{match});
    s = setfield(s,F{match},n);
  end
  if i==1, t=s; else n=s; end
end

catch me, excuse(me); end
