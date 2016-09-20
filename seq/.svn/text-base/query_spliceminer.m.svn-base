function res=query_spliceminer(gene)
% (Gaddy's code from March2008)
% http://discover.nci.nih.gov/spliceminer/Batch?queryBy=Gene&organism=9606&text=APC1

s=urlread('http://discover.nci.nih.gov/spliceminer/Batch','get',{'queryBy','Gene','organism','9606','text',gene});
res={};
ds=dlmsep(s,char(10));
for i=1:length(ds)
  if ~isempty(ds{i})
    res{end+1}=regexprep(dlmsep(ds{i},char(9)),'[\b\r\t ]','');
  end
end
if strcmp(res{2},'Noresultsfoundforthesymbolsintheuploadfile.')
  res={};
else
  res=cat(1,res{:});
end
