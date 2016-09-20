function b = parsehexcolors(a)

a = regexprep(a,'^#','');
a = regexprep(a,'^0x','');

b = nan(size(a,1),3);
for i=1:size(a,1)
  h = lower(a{i}); if length(h)~=6, error('hex colors should be 6-character strings'); end
  for j=1:length(h), z(j) = find(h(j)=='0123456789abcdef')-1; end
  b(i,:) = [z(1)*16+z(2) z(3)*16+z(4) z(5)*16+z(6)]/256;
end




