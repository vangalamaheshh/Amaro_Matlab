function p = read_params_file(param_file_name)
% Take in parameter file name, return a structure containing the parsed parameter values.

f=read_dlm_file(param_file_name);
for i=1:length(f)
  p(i).module=f{i}{1};
  p(i).param=f{i}{2};
  p(i).value=f{i}{3};
end
