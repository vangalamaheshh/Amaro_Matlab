function C=order_by_pos(C,secondary_sort)
% C = order_by_pos(C,secondary_sort)

if ~isfield(C,'chrn')
  C=add_chrn(C);
end

C.pos=as_column(C.pos);
C.chrn=as_column(C.chrn);

pos=double(C.chrn).*10.^(ceil(log10(max(double(C.pos))))+2)+double(C.pos);
if exist('secondary_sort','var') && ~isempty(secondary_sort)
  [sp,si]=sortrows([pos secondary_sort]);
else
  [sp,si]=sort(pos);
end

C=reorder_D_rows(C,si);
