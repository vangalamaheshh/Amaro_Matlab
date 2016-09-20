function m = listmap(a,b)

if ischar(a) && ischar(b) && length(b)<length(a)
  % special case
  method = 1;
else
  if ischar(a), a={a}; end
  if ischar(b), b={b}; end
  method = 2;
end

if method==3     % modest speedup with respect to method 2

 % STILL NOT WORKING

  m = nan(length(a),1);
  [tmp tmp aj] = unique(a);
  [tmp ia ib] = intersect(a,b);
  % first take care of all the singleton matches
  m(ia) = ib;
  % next take care of the non-unique matches
  h = histc(aj,1:length(a));
  idx = find(h>1);
  for j=1:length(idx), if ~mod(j,1e4), fprintf('%d/%d ',j,length(idx)); end
    i=idx(j);
    m(aj==i) = ib(i);
  end, if length(idx)>=1e4, fprintf('\n'); end

elseif method==2  % newer method, in use til Feb2013

  m = nan(length(a),1);
  [a ai aj] = unique(a);
  [c ia ib] = intersect(a,b);
  for i=1:length(c), if ~mod(i,1e5), fprintf('%d/%d ',i,length(c)); end
    m(aj==ia(i)) = ib(i);
  end, if length(c)>=1e5, fprintf('\n'); end

elseif method==1  % original method

  m = nan(length(a),1);
  [u ui uj] = unique(a);
  if ischar(u)
    for i=1:length(u)
      idx = find(b==u(i));
      if isempty(idx), idx = NaN; end
      m(uj==i) = idx;
    end
  else
    for i=1:length(u)
      idx = find(strcmp(u{i},b),1);
      if isempty(idx), idx = NaN; end
      m(uj==i) = idx;
    end
  end  
end
