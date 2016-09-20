function As = enum_mat_with_fixed_margins(R,C)
  % Function enumerates all m x n matrices with fixed row and column
  % marginal sums contained in the m x 1 vector R and the n x 1 vector C
    
  m = size(R,1);
  n = size(C,2);
  
  for i=1:length(R)
    Ps = enumerate_partitions(R(i));
    perm_R{i} = [];
    new_R = [];
    for j=1:n
      cur_Ps = find(cellfun(@length,Ps) == j);
      cur_R = [cat(2,Ps{cur_Ps})' zeros(length(cur_Ps),n-j)];
      new_R = [new_R; cur_R];
    end
    num_perms(i) = size(new_R,1);
    % now expand to include all unique 
    for k=1:num_perms(i)
      cur_perm = new_R(k,:);
      perm_R{i} = [perm_R{i}; unique(perms(cur_perm),'rows')];
    end
  end

  num_mats = prod(cellfun(@length,perm_R));
  
  
  row = mat2cell(perm_R{1},ones(size(perm_R{1},1),1),n);;
  for i=2:length(perm_R)
    row = merge_rows(row,mat2cell(perm_R{i},ones(size(perm_R{i},1),1),n));
  end
    
  % Now find only those matrices to keep -- those with column sums C    
  col_sums = cellfun(@sum,row,'UniformOutput',false);
  F = @(xx) prod(double(eq(xx,C)));
  
  keep = find(cellfun(F,col_sums));
  
  As = cell(1,length(keep));
  
  for j=1:length(keep)
    As{j} = row{keep(j)};
    
  end
  
function merged = merge_rows(row1,row2)
  
  merged = {};
  for i=1:length(row1)
    for j=1:length(row2)
      merged{length(merged)+1} = [row1{i}; row2{j}];
    end
  end
    
    
