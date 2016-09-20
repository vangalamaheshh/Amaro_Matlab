function ct = histcc(vals,bins,varargin)

ct = histc(vals,bins,varargin{:});

bintxt = num2cellstr(as_column(bins));

brack1 = repmat({'['},length(bins),1);
comma = repmat({','},length(bins),1);
brack2 = [repmat({')'},length(bins)-1,1);']'];

pr(['<';brack1;' '],[bintxt([1 1:end]);' '],[' ';comma;'>'],[' ';bintxt([2:end end end])],[' ';brack2;' '],...
 repmat({' '},length(bins)+2,1),[sum(vals<bins(1));as_column(ct);sum(vals>bins(end))]);

if nargout==0
  clear ct
end
