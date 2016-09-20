function out = expand_histogram(counts,values)

z = length(counts);
if z<numel(counts), error('matrix input not supported'); end

if exist('values','var')
  zz = length(values);
  if zz<numel(zz), error('matrix input not supported'); end
  if zz~=z, error('"values" should be same length as "counts"'); end
else
  values = 1:z;
end

out = cell(z,1);
for i=1:z
  out{i} = values(i)*ones(counts(i),1);
end

out = cat(1,out{:});

