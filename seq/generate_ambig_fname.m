function fname=generate_ambig_fname(i,freeze)
if ~exist('freeze','var') || isempty(freeze), freeze=2; end

lanes = load_lanelist;
freezes = load_freezelist;

if freeze<1 || freeze>slength(freezes), error('No such freeze'); end

fname = sprintf('/xchip/tcga/gbm/analysis/lawrence/wgs/unionambig/ambig_%03d.txt',i);
