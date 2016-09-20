function force_link(src,dest)
% force_link(src,dest)

fprintf('Linking %s -> %s\n',src,dest);
[a b] = system(['ln -sf ' src ' ' dest]);
if a>0
  error('Link failed');
  b
end
