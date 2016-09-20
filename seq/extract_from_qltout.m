function A = extract_from_qltout(file,ids)

in_file = file;
r = num2str(rand);
tmp_file = ['/xchip/tcga/gbm/analysis/lawrence/tmp/tmp_' r '.txt'];

problem_flag = false;
A = cell(length(ids),1);
for i=1:length(ids)
  id = ids(i);
  cmd = ['grep -P -m 1 "^QUERY\t' num2str(id) '\t" ' file ' > ' tmp_file];
  result = system(cmd);
  if result ~= 0
    fprintf('Problem extracting %d from %s\n',id,file);
    problem_flag = true;
    A{i} = [];
  else
    A{i} = parse_qltout(tmp_file);
    cmd = ['rm ' tmp_file];
    system(cmd);
  end
end
if ~problem_flag
  A = combine_structs(A);
else
  A = struct('foo',[],'bar',[]);
end
