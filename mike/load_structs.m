function x = load_structs(fname, varargin)

if ischar(fname)
  z = direc(fname);
elseif iscell(fname)
  z = fname;
else
  error('Unknown input argument type\n')
end

nz = length(z);
if nz==0, error('No such files(s): %s',fname); end

demand_files(z);

x = cell(length(z),1);
fprintf('Loading:\n');
for i=1:nz
  fprintf('\t%d/%d   %s\n',i,nz,z{i});
  x{i} = load_struct(z{i},varargin{:});
  x{i}.load_structs_filename = repmat({z{i}},slength(x{i}),1);
end
x = concat_structs_keep_all_fields(x);

x = order_fields_first(x,{'load_structs_filename'});


