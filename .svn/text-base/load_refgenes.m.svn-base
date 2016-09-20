function rg=load_refgenes(fname)

f=read_dlm_file(fname);
rg=[];
for i=2:length(f)
  rg(i-1).refseq=f{i}{1};
  rg(i-1).gene=f{i}{2};
  rg(i-1).locus_id=f{i}{3};
  rg(i-1).chr=f{i}{4};
  rg(i-1).strand=convert_enum(f{i}{5},{'+',1;'-',0});
  rg(i-1).start=str2num(f{i}{6});
  rg(i-1).end=str2num(f{i}{7});  
end

