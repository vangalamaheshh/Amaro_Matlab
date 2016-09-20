function S = tab2struct(filename, format, header_lines)
% alias for load_struct

if ~exist('header_lines', 'var')
  header_lines = 1;
end

if ~exist('lowercase_fieldnames', 'var')
  lowercase_fieldnames = false;
end

if ~exist('format', 'var')
  format = [];
end

S = load_struct(filename, format, header_lines);

