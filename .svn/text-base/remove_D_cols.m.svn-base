function [D,rm_idx]=remove_D_cols(D,match_str,varargin)

rm_idx=strmatch(match_str,D.sdesc,varargin{:});
D=reorder_D_cols(D,setdiff(1:size(D.dat,2),rm_idx));

