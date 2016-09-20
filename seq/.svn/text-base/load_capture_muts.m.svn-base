function M = load_capture_muts(direc)

d = dir([direc '/TCGA-*']);
M = {}; np=0;
for i=1:length(d)
  fname = [direc '/' d(i).name '/mutation_reports5f.maf.annotated'];
  if exist(fname,'file')
    np=np+1;
    M{np} = load_struct(fname);
  end
end
M = concat_structs(M(:));
