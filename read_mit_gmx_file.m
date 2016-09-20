function sets=read_mit_gmx_file(fname)

f=read_dlm_file(fname);

nsets=length(f{1});
for i=1:nsets
  sets{i}.name=f{1}{i};
  if isempty(sets{i}.name)
    nsets=i-1;
    disp(['Name of set no. ' num2str(i) ...
          ' is blank... stopping reading sets']);
    sets=sets(1:nsets);
    break;
  end
end

for j=1:nsets
  if ~strcmp(f{2}{j},'BLACK')
    disp(['Missing BLACK in set ' sets{j}.name ', column ' num2str(j)]);
  end
end

for i=3:length(f)
  for j=1:nsets
    st=f{i}{j};
    if ~strcmp(st,'null')
      sets{j}.genes{i-2}=st;
    end
  end
end
