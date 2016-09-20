function t = find_stretches(f)
  f = as_column(f);
  k = [0;find(diff(f))]+1;
  t = [f(k) k [k(2:end)-1;length(f)]];

  
