function [row filled] = rowperm(N,R,mx_row,row_ids)
  
  row = zeros(1,N);
  ids = 1:N;
  
  N1 = N - R(mx_row);
  
  cur_max = N;
  
  for i=1:length(row_ids)
    cur_row = row_ids(i);
    for j=1:R(cur_row)
      cur_id = ids(ceil(rand()*cur_max));
      row(cur_id) = cur_row;
      ids(cur_id) = cur_max;
      cur_max = cur_max-1;
    end
  end
  
  filled = find(row);
  row(find(row == 0)) = mx_row;
  
  
