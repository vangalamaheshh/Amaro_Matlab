function C=add_batch_id(C)

if ~iscell(C.sdesc)
    C.sdesc=cellstr(C.sdesc);
end

for i=1:length(C.sdesc)
  st{i}=regexprep(C.sdesc{i},'(.*)_','');
end

u=unique(strvcat(st),'rows');

e=convert_enum(deblank(st),[ cellstr(u) mat2cell((1:size(u,1))',ones(size(u,1),1),1) ]);

C=add_D_sup(C,'BATCH','Batch',e,'cols');
