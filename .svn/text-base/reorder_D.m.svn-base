function D=reorder_D(D,rowcol,idx)

switch rowcol(1:min(length(rowcol),3))
 case {'col','sam','con','exp'}
  D=reorder_D_cols(D,idx);
 case {'row','gen','mir'}
  D=reorder_D_rows(D,idx);
 otherwise
  error('no match');
end

