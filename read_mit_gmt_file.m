function sets=read_mit_gmt_file(fname)

fid=fopen(fname,'r');

ln=1;
while 1
  tline = fgetl(fid);
  if ~ischar(tline), break, end
  l{ln}=tline;
  ln=ln+1;
  if mod(ln,1000)==0
    disp(ln)
  end
end
fclose(fid);
ln=ln-1;

k=1;
for i=1:ln
  if isempty(l{i})
    disp(['Empty line ' num2str(i) ]);
    continue;
  end
  sd=dlmsep(l{i});
  cur_set=[];
  cur_set.name=sd{1};
  cur_set.type=sd{2};
%  if ~strcmp(sd{2},'BLACK')
%    disp(['Missing BLACK in set ' cur_set.name ]);
%  end
  for j=3:length(sd)
    cur_set.genes{j-2}=sd{j};
  end
  sets(k)=cur_set;
  k=k+1;
end



