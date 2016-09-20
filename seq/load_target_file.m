function T = load_target_file(fname,P)

if ~exist('P','var'), P=[]; end

demand_file(fname);

tmp = load_lines(fname);
numlines = length(tmp);
clear tmp;

fprintf('Loading target file: %s\n',fname);
cc = get_colcount(fname);
T = load_struct(fname,['%s %s' repmat('%f',1,cc-2)],char(9),0);

if strcmpi(T.col1{1},'gene')
  error('target list should LACK a header');
end

if slength(T)<numlines
  error('Problem loading target file!  Only %d/%d lines were able to be loaded.',slength(T),numlines);
end

if ~isfield(T,'col5'), error('Need GC content for capture plot!'); end
T = rename_field(T,colx(1:5),{'gene','chr','start','end','gc'});

T.chr = convert_chr(T.chr,P);

if isfield(T,'col6')
  T = rename_field(T,'col6','len');
else
  T.len = T.end-T.start+1;
end
if any(T.len ~= (T.end-T.start+1))
  error('"len" column in target file is calculated wrong');
end

if isfield(T,'col7')
  T = rename_field(T,'col7','membership');
end
if isfield(T,'col8')
  error('Extraneous column 8 in target file %s',fname);
end

