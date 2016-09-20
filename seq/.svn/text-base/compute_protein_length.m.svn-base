function len = compute_protein_length(code_start,code_end,exon_starts,exon_ends)
% Mike Lawrence 2010

if length(exon_starts)~=length(exon_ends), error('length(exon_starts)~=length(exon_ends)'); end

len = 0;
for i=1:length(exon_starts)
  s = exon_starts(i);
  e = exon_ends(i);
  if e<code_start || s>code_end, continue; end
  if s<code_start, s=code_start; end
  if e>code_end, e=code_end; end
  len = len + (e-s+1);
end

len = round(len/3);
