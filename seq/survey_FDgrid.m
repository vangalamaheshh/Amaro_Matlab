function FD = survey_FDgrid(dna,frame,nt5splice,nt3splice)
%
% Given a sequence of DNA, constructs an FDgrid, which is a matrix of the following form:
%
%                           col1   col2   col3   col4   col5   col6   col7   col8    col9
%                             A      T    C(CpG) C(TpC)   C    G(CpG) G(GpA)   G    unknown
%   row1 = synonymous         .      .      .      .      .      .      .      .      .
%   row2 = missense           .      .      .      .      .      .      .      .      .
%   row3 = nonsense           .      .      .      .      .      .      .      .      .
%   row4 = nonstop            .      .      .      .      .      .      .      .      .
%   row5 = splice-site        .      .      .      .      .      .      .      .      .
%   row6 = unknown            .      .      .      .      .      .      .      .      .
%
%   where the entire matrix sums to 3*length(dna), representing all possible mutations.
%
% "frame" = 0,1,or 2.
%    0 = first three bases are the first codon.
%    1 = skip 1 base, then the next three bases are the first codon.
%    2 = skip 2 bases, then the next three bases are the first codon.
%
% "nt5splice" and "nt3splice" define the number of 5' and 3' nucleotides to be counted
% as being at risk for splice-site mutations, a category which overrides all the others.
%
% col9 (unknown base composition) will always sum to (3 x two end positions) = 6
% row6 (unknown functional consequences) will sum to (3 x [at most 4]) depending on frame, ntXsplice
%
% to avoid nonzero unknowns, append flanking bases to the query sequence.
%
% Mike Lawrence 2008-08-04
%

if size(dna,1) ~= 1, error('"dna" must be a normal string'); end

if ~exist('frame','var'), frame=0; end
if ~exist('nt5splice','var'), nt5splice=0; end
if ~exist('nt3splice','var'), nt3splice=0; end

F = survey_all_mutations(dna,frame,nt5splice,nt3splice);
D = survey_contexts(dna);
D(D==0) = 9;     % unknown base context at edges

FD = zeros(6,9);
for i=1:length(dna)
  FD(:,D(i)) = FD(:,D(i)) + F(:,i);
end

