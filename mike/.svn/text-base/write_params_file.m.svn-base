function write_params_file(P,fname)

x=[];
x.key = fieldnames(P);
x.value = cell(length(x.key),1);
for i=1:length(x.key)
  z = getfield(P,x.key{i});
  if ischar(z), x.value{i} = z;
  elseif isnumeric(z)||islogical(z), x.value{i} = num2str(z);
  else 
    error('unknown parameter type for P.%s\n',x.key{i});
  end
end

save_struct_noheader(x,fname);
