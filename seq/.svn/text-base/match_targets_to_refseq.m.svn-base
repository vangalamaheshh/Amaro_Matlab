function match_targets_to_refseq(infile,R)
% match_targets_to_refseq(infile,R)
%
% Mike Lawrence 2010-02-17

if ~exist('R','var'), build = 'hg18', R = load_refseq(build); end
R.chrname = R.chr;
R.chr = convert_chr(R.chrname);
R = reorder_struct(R,~isnan(R.chr));

T = load_target_interval_list(infile)
T = reorder_struct(T,~isnan(T.chr));

% first convert refseq to exon table
% chr st en refseq_idx exon_no

e = cell(slength(R),1);
for i=1:slength(R)
  ne = R.n_exons(i);
  e{i} = [R.chr(i)*ones(ne,1) R.exon_starts{i} R.exon_ends{i} i*ones(ne,1) (1:ne)'];
end
e = cat(1,e{:});
e = sortrows(e);


%% (not finished)
