function gcc = gc_content(in)

if ~iscell(in), in = {in}; end

gcc = nan(length(in),1);
for i=1:length(in)
  d = upper(in{i});
  at = sum(d=='A' | d=='T');
  gc = sum(d=='G' | d=='C');
  gcc(i) = gc/(gc+at);
end
