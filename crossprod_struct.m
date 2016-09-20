function sa=crossprod_struct(sa1,sa2)

sa=add_struct(repmat(as_column(sa1),1,length(sa2)), ...
              repmat(as_row(sa2),length(sa1),1));
sa=sa(:);

