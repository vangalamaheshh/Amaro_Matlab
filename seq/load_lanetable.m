function T = load_lanetable(fname)

if ~exist(fname,'file'), error('%s not found',fname); end

try
  T = load_struct(fname,'%f%s%s%s%s%s%s%s');
catch me
  T = load_struct(fname,'%f%s%s%s%s%s%s');
  fprintf('No baitset column: filling in with "unknown"\n');
  T.baitset = repmat({'unknown'},slength(T),1);
end

tmp = parse(T.PU,'^(.....)....(\d\d\d\d\d\d)\.(\d)$',{'fcname','fcdate','fclane'});
T = merge_structs({T,tmp});
T = make_numeric(T,'fclane');
