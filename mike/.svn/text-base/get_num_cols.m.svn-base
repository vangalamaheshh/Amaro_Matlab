function n = get_num_cols(file)
% Mike Lawrence 2009-05-18

demand_file(file);

f = fopen(file);
nl=0;
for i=1:10
  tmp = fgetl(f);
  if tmp==-1, break; end
  l{i,1} = tmp;
  nl=nl+1;
end
fclose(f);
if nl==0, n=0; fprintf('Warning: empty file %s\n',file); return; end

for i=1:nl, nt(i,1) = sum(l{i}==char(9))+1; end

if all(nt==nt(1))
  n = nt(1);
else
  fprintf('WARNING: Inconsistent column counts in %s:\n',file);
  for i=1:nl, fprintf('  Line %d: %d columns\n',i,nt(i)); end
  n = max(nt);
  fprintf('Assuming file width = maximum of first ten lines = %d\n',n);
end

