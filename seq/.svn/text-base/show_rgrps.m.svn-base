function show_rgrps(f)
dsi = f.getFileHeader().getReadGroups().iterator();;
while (dsi.hasNext())
  s = dsi.next();
  name = s.getAttribute('PU'); % already char
  index = s.getReadGroupId().toCharArray();
  fprintf('%s\t%s\n',index,name);
end
