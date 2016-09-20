function write_graph_xml(fname,nodes,C,is_directed)
% write_graph
if ~exist('is_directed','var')
  is_directed=0;
end

f=fopen(fname,'w');
fprintf(f,'<!-- Matlab graph writer :: %s -->\n',datestr(now));
fprintf(f,'<graph directed="%d">\n',is_directed);

fprintf(f,'   <!-- nodes -->\n');
x=500+400*sin((0:(length(nodes)-1))/length(nodes)*2*pi);
y=500+400*cos((0:(length(nodes)-1))/length(nodes)*2*pi);
for i=1:length(nodes)
  fprintf(f,'   <node id="%s" label="%s">\n',nodes(i).id,nodes(i).label);
  fprintf(f,'      <att name="Y" value="%5.1f"/>\n',y(i));
  fprintf(f,'      <att name="X" value="%5.1f"/>\n',x(i));
  fprintf(f,'   </node>\n');
  % add additional attributes
end

fprintf(f,'   <!-- edges -->\n');
if ~is_directed
  C=triu(C,1);
end
[source,dest,value]=find(C);

for i=1:length(source)
  fprintf(f,'   <edge source="%s" target="%s">\n',nodes(source(i)).id,nodes(dest(i)).id);
  fprintf(f,'      <att name="weight" value="%f"/>\n',abs(value(i)));
  fprintf(f,'      <att name="sign" value="%d"/>\n',sign(value(i)));
  fprintf(f,'   </edge>\n');
end
fprintf(f,'</graph>\n');
fclose(f);
