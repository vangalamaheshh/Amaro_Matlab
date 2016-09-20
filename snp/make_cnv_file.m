function make_cnv_file(fname,outfname,GI_fname,columns)

if ~exist('columns','var')
  columns=[2 3 4 5];
end
disp([ 'Using columns: ' num2str(columns)]);


GI_fname=regexprep(GI_fname,'\.txt$','');
if exist([GI_fname '.mat'])
  load([GI_fname '.mat']);
else
  GI=read_genome_info_file([ GI_fname '.txt']);
  GI.dat=cat(1,GI.dat{:});
end
if nnz(cellfun('isempty',GI.dat(:,2)))>0
  error('has empty cells');
end
chrn=chromosome2num(GI.dat(:,2));

if nnz(cellfun('isempty',GI.dat(:,4)))>0
  error('has empty cells');
end
pos=str2num(strvcat(GI.dat(:,4)));

d=read_dlm_file(fname);
t={};
for i=2:length(d)
  id=d{i}{columns(1)};
  c=chromosome2num(d{i}{columns(2)});
  st=str2num(d{i}{columns(3)});
  en=str2num(d{i}{columns(4)});
  idx=find(chrn==c & pos>=st & pos<=en);
  if ~isempty(idx)
    t=[t; GI.dat(idx,1) cellstr(repmat(id,length(idx),1))];
  end
end

write_table(outfname,t);
