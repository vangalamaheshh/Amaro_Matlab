function T = calculate_TCC_matrix(ref,targ)
% calculate_TCC_matrix(ref,targ)
%
% inputs:
%
%   ref = reference transcript model
%      with fields: build, chr, strand, coding_start, coding_end, exon_starts, exon_ends
%   targ = target region list
%      with fields: build, chr, region_starts, region_ends
%
% outputs:
%
%   12x16x6 TCC matrix =
%
%     12 rows for mutation Type (e.g. A->C)
%        A->C = 1      A->G = 2      A->T = 3
%        C->A = 4      C->G = 5      C->T = 6
%        G->A = 7      G->C = 8      G->T = 9
%        T->A = 10     T->C = 11     T->G = 12
%
%     16 columns for DNA Context (e.g. A_C) w.r.t. (+) strand
%        A_A = 1     A_C = 2     A_G = 3     A_T = 4
%        C_A = 5     C_C = 6     C_G = 7     C_T = 8
%        G_A = 9     G_C = 10    G_G = 11    G_T = 12
%        T_A = 13    T_C = 14    T_G = 15    T_T = 16
%
%     6 pages for Consequence
%        1 = synonymous
%        2 = missense
%        3 = nonsense
%        4 = nonstop
%        5 = splice_site (first & second bases outside exon)
%        6 = UTR / intron / IGR / unknown
%
% Mike Lawrence 2009-01-21

% STEP 0:  Check input

require_fields(ref,{'build','chr','strand','coding_start','coding_end','exon_starts','exon_ends'});
require_fields(targ,{'build','chr','region_starts','region_ends'});

if ~strcmp(ref.build,targ.build), error('Builds do not agree'); end

if isnumeric(ref.chr), ref.chr = num2str(ref.chr); end;
if isnumeric(targ.chr), ref.targ = num2str(targ.chr); end;
if ~strncmp(ref.chr,'chr',3), ref.chr = ['chr' ref.chr]; end;
if ~strncmp(targ.chr,'chr',3), targ.chr = ['chr' targ.chr]; end;
if ~strcmp(ref.chr,targ.chr), error('Chromosomes do not agree'); end

ref.start = min(ref.exon_starts);
ref.end = max(ref.exon_ends);
targ.start = min(targ.region_starts);
targ.end = max(targ.region_ends);
if ref.end < targ.start || ref.start > targ.end, error('Target and reference do not overlap'); end

% STEP 1:  Retrieve DNA sequence of entire target

dna_min = min(ref.start,targ.start)-2;
dna_max = max(ref.end,targ.end)+2;
dna = genome_region(ref.chr,dna_min,dna_max,ref.build);
dna_offset = dna_min-1;
dna = upper(dna);

% STEP 2:  Transcribe to mRNA sequence

rna = [];
rna2dna = [];

for i=1:length(ref.exon_starts)
  st = ref.exon_starts(i);
  en = ref.exon_ends(i);
  if ref.coding_start<en && ref.coding_end>st
    st = max(st,ref.coding_start);
    en = min(en,ref.coding_end);
    rna = [rna dna(st-dna_offset:en-dna_offset)];
    rna2dna = [rna2dna st:en];
  end
end

minusstrand = strcmp(ref.strand,'-');
if minusstrand
  rna = my_seqrcomplement(rna);
  rna2dna = fliplr(rna2dna);
end

% truncate incomplete terminal codon

x = 3*floor(length(rna)/3);
rna = rna(1:x);
rna2dna = rna2dna(1:x);

% mark correspondence to DNA

dna2rna = nan(1,length(dna));
dna2rna(rna2dna-dna_offset) = 1:length(rna);

% STEP 3:  Calculate consequences of mutations in mRNA

csq = zeros(3,length(rna));

cpos = 0;
bases = 'ACGT';
for i=1:length(rna)
  cpos=cpos+1;
  if cpos==4, cpos=1; end
  if cpos==1, old_codon = rna(i:i+2); old_aa = my_nt2aa(old_codon); end
  new_bases = setdiff(bases,rna(i));
  for j=1:3
    new_codon = old_codon;
    new_codon(cpos) = new_bases(j);
    new_aa = my_nt2aa(new_codon);
    if strcmp(old_aa,new_aa), csq(j,i) = 1;  % synonymous
    elseif strcmp(new_aa,'*'), csq(j,i) = 3; % nonsense
    elseif strcmp(old_aa,'*'), csq(j,i) = 4; % nonstop
    else csq(j,i) = 2; end                   % missense
  end
end

% STEP 4:  Extract the target positions
%    For each position in the target region
%        Determine Context
%        For each of the three changes
%          Determine Type
%          If it's in the RNA, read Consequence from csq
%          Else choose 5 or 6 based on proximity to splice-site

% compute non-overlapping targets
tdna = zeros(dna_max-dna_offset,1);
for t=1:length(targ.region_starts)
  tdna(targ.region_starts(t)-dna_offset:targ.region_ends(t)-dna_offset) = 1;
end

T = zeros(12,16,6);
for i=dna_min:dna_max
  if tdna(i-dna_offset)
    if ~minusstrand
      flank5 = dna((i-1)-dna_offset);
      old_base = dna((i)-dna_offset);
      flank3 = dna((i+1)-dna_offset);
    else
      flank5 = my_seqrcomplement(dna((i+1)-dna_offset));
      old_base = my_seqrcomplement(dna((i)-dna_offset));
      flank3 = my_seqrcomplement(dna((i-1)-dna_offset));
    end
    flank5 = find(bases==flank5)-1;
    old_base = find(bases==old_base)-1;
    flank3 = find(bases==flank3)-1;
    context = (flank5*4)+1 + flank3;
    rna_pos = dna2rna(i-dna_offset);
    if ~isnan(rna_pos)
      consequences = csq(:,rna_pos);
    else
      consequences = [6;6;6];      % unknown/UTR/intron/IGR
      for ii=[i-2 i-1 i+1 i+2]
        if ~isnan(dna2rna(ii-dna_offset))
          consequences = [5;5;5];  % splice-site
          break;
    end,end,end
    for j=1:3
      type = (old_base*3)+j;
      row = type;
      col = context;
      page = consequences(j);
      T(row,col,page) = T(row,col,page) + 1;

if 0
  targ
  [targ.region_starts targ.region_ends]
  ref
  [ref.exon_starts ref.exon_ends] 
  i
  T(row,col,page) = -T(row,col,page); 
  T
  T(row,col,page) = -T(row,col,page);
  keyboard
end

    end
  end
end
      

