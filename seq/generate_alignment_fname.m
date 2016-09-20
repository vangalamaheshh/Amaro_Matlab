function fname=generate_alignment_fname(i,freeze)
if ~exist('freeze','var') || isempty(freeze), freeze=2; end
fname = generate_qltout_fname(i,freeze);
