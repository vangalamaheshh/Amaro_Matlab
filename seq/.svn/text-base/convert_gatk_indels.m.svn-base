function M = convert_gatk_indels(M)
% convert_gatk_indels(M)
%  demand_fields(M,{'start','end','ref_allele','newbase'});
%
% DELETIONS:
%    finds all records where length(newbase) < length(ref_allele) 
%       -- removes redundant bases and adjusts start,end
%              e.g. "TC to T" changed to "C to -" and start++, end++
%                   "TTGA to "T" changed to "TGA to -" and start++, end++
%                   "CTGA to "CT" changed to "GA to -" and start+=2, end+=2
%
% INSERTIONS:
%    finds all records where length(newbase) > length(ref_allele)
%       -- removes redundant bases and adjusts start,end
%              e.g. "T to TC" changed to "- to C", no change to start/end
%                   "CGTAAA" to "CGTAAAC" changed to "- to C" and start+=5, end+=5
%                   "A to AGGG" changed to "- to GGG", no change to start/end

demand_fields(M,{'start','end','ref_allele','newbase'});
M = make_numeric(M,{'start','end'});
M.ref_allele = upper(M.ref_allele);
M.newbase = upper(M.newbase);

M.reflen = cellfun('length',M.ref_allele);
M.reflen(strcmp(M.ref_allele,'-')) = 0;
M.altlen = cellfun('length',M.newbase);
M.altlen(strcmp(M.newbase,'-')) = 0;

ins = find(M.altlen>M.reflen & M.reflen>0);
del = find(M.altlen<M.reflen & M.altlen>0);

%DELETIONS
for j=1:length(del), i=del(j);
  ref = M.ref_allele{i};
  alt = M.newbase{i};
  a = length(alt);
  z = find(ref(1:a)~=alt);
  if isempty(z)
    M.ref_allele{i} = M.ref_allele{i}(a+1:end);
    M.newbase{i} = '-';
    M.start(i) = M.start(i) + a;
    M.end(i) = M.end(i) + a;
  else
    fprintf('Complex deletion left as is:  %s to %s\n', ref,alt); 
  end
end

%INSERTIONS
for j=1:length(ins), i=ins(j);
  ref = M.ref_allele{i};
  alt = M.newbase{i};
  a = length(ref);
  z = find(alt(1:a)~=ref);
  if isempty(z)
    M.newbase{i} = M.newbase{i}(a+1:end);
    M.ref_allele{i} = '-';
    M.start(i) = M.start(i) + (a-1);
    M.end(i) = M.end(i) + (a-1);
  else
    fprintf('Complex insertion left as is:  %s to %s\n', ref,alt);
  end
end
