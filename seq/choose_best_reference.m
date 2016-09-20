function r = choose_best_reference(refs,targ)
% choose_best_reference(refs,targ)
%
% inputs:
%
%   refs = cell-array of reference transcript models
%      each cell having fields: build, chr, strand, coding_start, coding_end, exon_starts, exon_ends
%   targ = target region list
%      having fields: build, chr, region_starts, region_ends
%
% output:
%
%   r = the index of the reference that best fits the target region list
%
% procedure:
%
%   for each reference transcript model, compares to target list and computes the following:
%     to = number of target nucleotides missing from the reference(coding sequence)
%     ro = number of reference(coding sequence) nucleotides missing from the target
%     penalty = 4*to+ro
%   returns the r with lowest penalty
%
% Mike Lawrence 2009-01-22

require_fields(targ,{'build','chr','region_starts','region_ends'});
targ.start = min(targ.region_starts);
targ.end = max(targ.region_ends);

if ~iscell(refs), error('refs should be a cell array'); end
nr = length(refs);

to = nan(nr,1);
ro = nan(nr,1);

for i=1:nr
  require_fields(refs{i},{'build','chr','strand','coding_start','coding_end','exon_starts','exon_ends'});
  if ~strcmp(refs{i}.build,targ.build), error('Builds do not agree for ref %d',i); end
  if isnumeric(refs{i}.chr), refs{i}.chr = num2str(refs{i}.chr); end;
  if isnumeric(targ.chr), refs{i}.targ = num2str(targ.chr); end;
  if ~strncmp(refs{i}.chr,'chr',3), refs{i}.chr = ['chr' refs{i}.chr]; end;
  if ~strncmp(targ.chr,'chr',3), targ.chr = ['chr' targ.chr]; end;
  if ~strcmp(refs{i}.chr,targ.chr), error('Chromosomes do not agree for ref %d',i); end

  refs{i}.start = min(refs{i}.exon_starts);
  refs{i}.end = max(refs{i}.exon_ends);

  dna_min = min(refs{i}.start,targ.start)-2;
  dna_max = max(refs{i}.end,targ.end)+2;
  dna_offset = dna_min-1;
  d = zeros(dna_max-dna_offset,1);
  for j=1:length(targ.region_starts)
    d(targ.region_starts(j)-dna_offset:targ.region_ends(j)-dna_offset) = 1;
  end
  for j=1:length(refs{i}.exon_starts)
    st = refs{i}.exon_starts(j);
    en = refs{i}.exon_ends(j);
    if refs{i}.coding_start<en && refs{i}.coding_end>st
      st = max(st,refs{i}.coding_start);
      en = min(en,refs{i}.coding_end);
      d(st-dna_offset:en-dna_offset) = d(st-dna_offset:en-dna_offset) + 2;
    end
  end
  to(i) = sum(d==1);
  ro(i) = sum(d==2);
end

penalty = 4*ro + to;
[tmp,r] = min(penalty);
