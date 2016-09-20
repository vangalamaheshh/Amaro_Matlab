function X = load_struct_specify_numeric_cols(fname,numeric_cols,num_header_lines,lowercase_fieldnames)
% load_struct_specify_numeric_cols(fname,numeric_cols,num_header_lines)
%
% Mike Lawrence 2009-07-01

if ~exist('fname','var'), error('Must specify fname'); end
if ~exist('numeric_cols','var'), numeric_cols = []; end
if ~exist('num_header_lines','var'), num_header_lines = 1; end
if ~exist('lowercase_fieldnames','var'), lowercase_fieldnames = false; end

if ~exist(fname,'file'), error('%s not found',fname); end

f = fopen(fname);
for i=1:num_header_lines+1; l = fgetl(f); end
numcols = sum(l==char(9))+1;

is_string = true(1,numcols);
is_string(numeric_cols) = false;
format = [];
for i=1:numcols
  if is_string(i), format = [format '%s'];
  else format = [format '%f']; end
end

%X = load_struct(fname,format,num_header_lines,lowercase_fieldnames);
P=[]; P.lowercase_fieldnames = lowercase_fieldnames;
X = load_struct(fname,format,char(9),num_header_lines,P);

