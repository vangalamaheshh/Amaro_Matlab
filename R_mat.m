function r=R_mat(x)

r='rbind(';
for i=1:size(x,1)
  s=num2str(x(i,:),'%d,');
  r=[r 'c(' s(1:(end-1)) '),'];
end
r=[r(1:(end-1)) ')'];

