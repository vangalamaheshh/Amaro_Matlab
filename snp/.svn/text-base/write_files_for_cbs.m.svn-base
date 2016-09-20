function write_files_for_cbs(X,dirname,prefix)

warning(['This function has been deprecated. Please use run_cbs ' ...
         'instead. NS 2009.02.09'] )

cd(dirname);

f=cell(size(X.dat,2),1);
fname=cell(size(X.dat,2),1);

for i=1:size(X.dat,2)
  fname{i}=[ prefix '_Sample' num2str(i)];
  f{i}=fopen(fname{i},'w');
  fprintf(f{i},'Marker\tChromosome\tPhysicalPosition\t%s\n',X.sdesc{i});
end

disp('taking 2.^(X+1)');
X.dat=2.^(X.dat+1);
disp('Writing files');
if ~isfield(X,'chr')
  X.chr=num2chromosome(X.chrn);
end
chunk=50000;
for i=1:chunk:size(X.dat,1)
  cur_idx=i:min(i+chunk-1,size(X.dat,1));
  c_marker=as_row(X.marker(cur_idx));
  c_chr=as_row(X.chr(cur_idx));  
  c_pos=mat2cell(as_row(X.pos(cur_idx)),1,ones(length(cur_idx),1));
  c=[c_marker; c_chr; c_pos];
  for j=1:size(X.dat,2)
    c_dat=mat2cell(X.dat(cur_idx,j)',1,ones(length(cur_idx),1));
    cc=[c; c_dat];
    fprintf(f{j},'%s\t%s\t%d\t%0.3f\n',cc{:});
  end
  fprintf(1,'.%d',max(cur_idx));
end

for i=1:size(X.dat,2)
  fclose(f{i});
end
