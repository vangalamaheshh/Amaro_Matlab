function dRanger_write_links(fname,t)

x = [];
n=1;
for i=1:length(t)
  for j=1:size(t{i}.annotated_support,1)
    x = [x sprintf('link%05d hs%d %d %d\nlink%05d hs%d %d %d\n',...
       n,t{i}.annotated_support(j,7),t{i}.annotated_support(j,9),t{i}.annotated_support(j,10),...
       n,t{i}.annotated_support(j,11),t{i}.annotated_support(j,13),t{i}.annotated_support(j,14))];
    n = n + 1;
  end
end
save_textfile(x,fname);
