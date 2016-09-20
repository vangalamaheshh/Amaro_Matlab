function bainame = find_bai(bamname)

if ~iscell(bamname)
  if ~exist(bamname), error('bamfile %s does not exist', bamname); end
  bainame = regexprep(bamname,'\.bam$','\.bai');
  if strcmp(bainame,bamname) || ~exist(bainame,'file')
    bainame = [bamname '.bai'];
    if ~exist(bainame,'file'), error('bamfile needs to have an index in the same directory'); end
  end
else
  for i=1:length(bamname)
    if ~exist(bamname{i}), error('bamfile %s does not exist', bamname{i}); end
    bainame{i,1} = regexprep(bamname{i},'\.bam$','\.bai');
    if strcmp(bainame{i},bamname{i}) || ~exist(bainame{i},'file')
      bainame{i} = [bamname{i} '.bai'];
      if ~exist(bainame{i},'file'), error('bamfile needs to have an index in the same directory'); end
    end
  end
end
