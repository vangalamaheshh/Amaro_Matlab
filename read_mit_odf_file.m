function T=read_mit_odf_file(fname)
% T.n_header
% T.col_names
% T.col_types
% T.xxx
% T.data

T.rem={};
fid=fopen(fname,'r');
tline=fgetl(fid); % title name
tline='';
while isempty(tline) || tline(1)=='#'
  if ~isempty(tline) && tline(1)=='#'
    T.rem{end+1}=deblank(tline);
  end
  tline=fgetl(fid);
end
t=dlmsep(tline,'=');
T.n_header=sscanf(t{2},'%d');

tline='';
while isempty(tline) || tline(1)=='#' 
  if ~isempty(tline) && tline(1)=='#'
    T.rem{end+1}=deblank(tline);
  end
  tline=fgetl(fid); %COLUMN NAMES:
end

t=dlmsep(tline,':');
t=dlmsep(t{2},9);
if isempty(t{1})
  T.col_names=t(2:end);
else
  T.col_names=t;
end

tline='';
while isempty(tline) || tline(1)=='#' 
  if ~isempty(tline) && tline(1)=='#'
    T.rem{end+1}=deblank(tline);
  end
  tline=fgetl(fid); %COLUMN TYPES:
end
t=dlmsep(tline,':');
if ~isempty(t) && strcmp(t{1},'COLUMN_TYPES')
  t=dlmsep(t{2},9);
  if isempty(t{1})
    T.col_types=t(2:end);
  else
    T.col_types=t;
  end
  use_tline=0;
else
  T.col_types=cellstr(repmat('string',length(T.col_names),1))';
  use_tline=1;
end

for i=(3-use_tline):T.n_header
  if ~use_tline
    tline='';
    while isempty(tline) || tline(1)=='#' 
      if ~isempty(tline) && tline(1)=='#'
        T.rem{end+1}=deblank(tline);
      end
      tline=fgetl(fid);
    end
  end
  use_tline=0;
  t=dlmsep(tline,'=');
  T=setfield(T,fix_field_name(t{1}),t{2}); %sscanf(t{2},'%s'));
end

if isfield(T,'DataLines')
  T=rmfield(T,'DataLines');
end

form='';
bools=[];
for i=1:length(T.col_types)
  switch (lower(T.col_types{i}))
   case {'int','integer'}
    form=[form '%d'];
   case 'string'
    form=[form '%s'];
   case 'float'
    form=[form '%f'];
   case 'boolean'
    form=[form '%s'];
    bools=[bools i];
  end
end

T.data=textscan(fid,form,'treatAsEmpty','NA');
for i=1:length(bools)
    bi=bools(i);
    x=zeros(size(T.data{bi},1),1);
    x(strmatch('true',T.data{bi}))=1;
    T.data{bi}=x;
end 

fclose(fid);
