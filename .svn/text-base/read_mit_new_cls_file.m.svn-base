function D=read_mit_new_cls_file(D,fname)

f=read_dlm_file(fname);

% handle # at 1st char
remove_list=[];
for i=1:length(f)
  if isempty(f{i}) || isempty(f{i}{1}) || f{i}{1}(1)=='#'
    remove_list=[ remove_list i];
  end
end
f=f(setdiff(1:length(f),remove_list));
verbose(['removing commented lines: ' num2str(remove_list) ],10);

% make sure all rows have the same number of columns
try
  t=cat(1,f{:});
catch
  n=cellfun('length',f);
  disp(n);
  error('read_mit_new_cls_file',...
        'not all lines have the same number of fields');
end

% is there a type row?
if (size(t,1)>1) && strcmp(lower(t{2,1}),'type')
  coltype=t(2,2:end);
  has_type_line=1;
  t=t([1 3:end],:); % remove type line
else
  coltype=(cellstr(repmat('categoric',size(t,2)-1,1)))';
  has_type_line=0;
  verbose('assumes all columns are categoric',10);
end

% handle annotation titles
supdesc=deblank((t(1,2:end))');
supacc=supdesc;
t=t(2:end,:); 

% get sample names and compare to ones in D
sdesc=t(:,1);
t=t(:,2:end); % remove sample name col
[M,mi,mj]=match_string_sets(sdesc,D.sdesc);

if iscell(D.sdesc)
  nDs=length(D.sdesc);
else
  nDs=size(D.sdesc,1);
end

if length(unique(mj))~=nDs
  error('read_mit_new_cls_file',...
        'sample names dont match');
end

% reorder t to match D (if needed)
if range(mi-(1:length(mi))')==0
  verbose('same order as D - no need to reorder.');
else
  t=t(mi,:);
  sdesc=sdesc(mi);
  verbose('reordering to match D.');
end

% remove comment cols
comment_cols=strmatch('comment',lower(coltype),'exact');
if ~isempty(comment_cols)
  com_t=t(:,comment_cols);
  com_supacc=supacc(comment_cols);
  com_supdesc=supdesc(comment_cols);
  keep_cols=setdiff(1:size(t,2),comment_cols);
  coltype=coltype(keep_cols);
  supdesc=supdesc(keep_cols);
  supacc=supacc(keep_cols);
  t=t(:,setdiff(1:size(t,2),comment_cols));
  verbose(['removing comment column(s):' num2str(as_row(comment_cols))]); 
end

% replace null/na with NaN
na1=strmatch('null',lower(t));
na2=strmatch('na',lower(t));
na12=[na1; na2];
if ~isempty(na12)
  t(na12)=mat2cell(repmat('NaN',length(na12),1),ones(length(na12),1),3);
end

%fill supdat
supdat=zeros(size(t,2),size(t,1));
for i=1:size(t,2)
  switch lower(coltype{i})
   case 'categoric'    
    [u,ui,uj]=unique_keepord(strvcat(t(:,i)),'rows');
    supdat(i,:)=uj;
    supdesc{i}=[ supdesc{i} ': '];
    for j=1:length(ui)
      supdesc{i}=[ supdesc{i} num2str(j) '-' deblank(u(j,:)) ];
      if j<length(ui)
        supdesc{i}=[ supdesc{i} '/'];
      end
    end
    supacc{i}=supdesc{i};
   case 'numeric'
    supdat(i,:)=str2num(strvcat(t(:,i)));
  end
end
D.supdat=supdat;
D.supacc=strvcat(supacc);
D.supdesc=strvcat(supdesc);
if exist('com_t','var')
  D.comacc=com_supacc;
  D.comdesc=com_supdesc;
  D.comdat=com_t';
end
