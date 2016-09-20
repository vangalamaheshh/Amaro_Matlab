function gc = get_gc_content(chrs,starts,ends,build)

if length(chrs)~=length(starts) || length(chrs)~=length(ends), error('input length mismatch'); end
if ~isnumeric(chrs) || ~isnumeric(starts) || ~isnumeric(ends), error('inputs should be numeric'); end

if ~exist('build','var')
  error('Must supply build!');
%  fprintf('Assuming build = hg18\n');
%  build = 'hg18';
end

gc = nan(length(chrs),1);
for c=1:get_chrcount(build)
  idx = find(chrs==c);
  if isempty(idx), continue; end
  fprintf('chr%d\n',c);
  ref = upper(genome_region(c,1,inf,build));
  len = length(ref);
  gc2 = nan(len,1);
  gc2(ref=='A' | ref=='T') = 0;
  gc2(ref=='C' | ref=='G') = 1;
  for j=1:length(idx), i=idx(j);
    st = max(1,starts(i));
    en = min(len,ends(i));
    gc(i) = nanmean(gc2(st:en));
end,end, fprintf('\n');
