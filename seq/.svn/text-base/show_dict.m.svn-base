function show_dict(f)
dsi = f.getFileHeader().getSequences().iterator();
while (dsi.hasNext())
  s = dsi.next();
  name = s.getSequenceName().toCharArray();
  index = s.getSequenceIndex();
  fprintf('%d\t%s\n',index,name);
end
