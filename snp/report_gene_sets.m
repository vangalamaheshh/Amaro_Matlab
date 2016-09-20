function report_gene_sets(fname,D,set_ids,comp_ids)

f=fopen(fname,'w');
[d,res]=dist(D.gsupdat(comp_ids,:),D.gsupdat(set_ids,:),struct('method','fisher','fdr',0.05));

ad=abs(d);
for i=1:size(d,1)
  [ds,di]=sort(ad(i,:));
  nonz=find(ds>0);
  di=di(nonz);
  di=fliplr(di);
  fprintf(f,'%s\t%d\n',deblank(D.gsupdesc(comp_ids(i),:)),full(sum(D.gsupdat(comp_ids(i),:),2)));
  fprintf(f,'%s\n',repmat('-',length(deblank(D.gsupdesc(comp_ids(i),:))),1));
  for j=1:length(nonz)
    fprintf(f,'%d\t%s\t%d\t%f\t%d\t%d\t%d\t%d\n',set_ids(di(j)),deblank(D.gsupdesc(set_ids(di(j)),:)),...
            full(sum(D.gsupdat(set_ids(di(j)),:),2)),d(i,di(j)),res{1}(i,di(j)),res{2}(i,di(j)),res{3}(i,di(j)),res{4}(i,di(j)));
  end
  fprintf(f,'\n');
end

fclose(f);

