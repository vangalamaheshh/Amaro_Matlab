function s=reorder_struct_inv(s,order)
s = reorder_struct(s,setdiff(1:slength(s),order));
