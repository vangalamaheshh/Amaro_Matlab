function s = sa2s(sa)
% sa2s(sa)
%
% converts struct array to struct
%
f = fields(sa);
c = squeeze(struct2cell(sa))';
s = [];
for i=1:length(f)
  x = c(:,i);
  % try to convert cell array of numbers to matrix array of numbers
  y = zeros(size(x));
  fail = false;
  for j=1:length(x)
    if ~isnumeric(x{j}), fail = true; break; end;
    y(j) = x{j};
  end
  if ~fail, x=y; end
  s = setfield(s,f{i},x);
end
