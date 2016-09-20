function M = crosstab2d(col1,col2)
  
  if size(col1,1) ~= size(col2,1)
    error('columns must be same length');
  end
  
  m = max(col1);
  n = max(col2);
  
  idx = n*(col1-1)+col2;
  
  hc = histc(idx,1:max(idx));
    
  if length(hc) ~= m*n
    hc(length(hc)+1:m*n) = zeros(1,m*n-length(hc));
  end
  
  xx = mat2cell(hc',1,repmat(n,1,m));
  M = cat(1,xx{:});
  
