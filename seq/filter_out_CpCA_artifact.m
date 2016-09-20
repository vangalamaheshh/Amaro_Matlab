function [A] = filter_out_CpCA_artifact(A, build, d) 

if ~isstruct(A) 
  A = load_struct(A)
elseif ~ischar(A) & ~isstruct(A) 
  error('MAF file must be either a string as path to the file or a struct!\n')
end

CONTEXT_NUMS_PLUS = [21,22,23,24, 31];
CONTEXT_NUMS_MINUS = [35,39,43,47,37 ]; 
AF_CUTOFF = 0.1;

if ~isfield(A, 'Chromosome'), A.Chromosome = A.chr; end 
if ~isfield(A, 'Start_position'), A.Start_position = A.start; end
if ~isfield(A, 'newbase'), A.newbase = A.Tumor_Seq_Allele2; end


categdir = fullfile('/xchip/cga1/lawrence/db', build, 'context65');
A.Chromosome = convert_chr(A.Chromosome); 
A = make_numeric(A, {'Start_position', 'i_tumor_f'});
A.context65 = get_context(A.Chromosome, A.Start_position, categdir);


% (+) strand
%categno = 21;     % C*CA -> A
newbase = 'A';
idx = arrayfun(@(x) ismember(x, CONTEXT_NUMS_PLUS), A.context65);
idx4 = idx & strcmp(A.newbase,newbase) & A.i_tumor_f < 0.1;

%% (-) strand
%categno = 47;     % T*GG -> T
newbase = 'T';
idx = arrayfun(@(x) ismember(x, CONTEXT_NUMS_MINUS), A.context65);
idx5 = idx & strcmp(A.newbase,newbase) & A.i_tumor_f < 0.1;

% filter out 
if d ==1 
  A = reorder_struct(A, ~(idx4 | idx5));
elseif d == -1 
  A = reorder_struct(A, (idx4 | idx5));
end