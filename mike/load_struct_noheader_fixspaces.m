function X = load_struct_noheader(fname,numcols,colnames)
% load_struct_noheader(fname[,numcols][,colnames])
% or 
% load_struct_noheader(fname[,format][,colnames])

% Mike Lawrence 2010-01-20

if ~exist('fname','var'), error('Must specify fname'); end
demand_file(fname);

% first fix spaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = fopen(fname); l = fgetl(f); fclose(f);
if ~contains(l,char(9)) && contains(l,char(32))
  t = load_textfile(fname);
  t = regexprep(t,'(\S)\s+(\S)','$1\t$2');
  save_textfile(t,fname);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==2 && iscell(numcols)
  colnames = numcols;
  numcols = [];
end

if ~exist('numcols','var') || isempty(numcols)
  f = fopen(fname);
  l = fgetl(f);
  numcols = sum(l==char(9))+1;
end

if isnumeric(numcols)
  format = repmat('%s',1,numcols);
else
  format = numcols;
end
X = load_struct(fname,format,char(9),0);

if exist('colnames','var')
  nf = length(fieldnames(X));
  nc = length(colnames);
  n = min([nf,nc]);
  X = rename_fields(X,colx(1:n),colnames(1:n));
  X = orderfields_first(X,colnames(1:n));
end
