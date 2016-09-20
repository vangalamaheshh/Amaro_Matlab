function [D1m,D2m]=match_Ds(D1,D2)

[Mt,m1,m2]=match_string_sets(D1.sdesc,D2.sdesc);
if length(unique(m2))~=length(m2)
  error('non unique sdesc');
end

D1m=reorder_D_cols(D1,m1);
D2m=reorder_D_cols(D2,m2);

