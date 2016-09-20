function [si,supdesc,supacc,supdat]=read_dchip_info_file(fname)

f=read_dlm_file(fname);
colnames=f{1};
f=f(2:end);
n=length(colnames);
ec=find(cellfun('isempty',colnames));

r=setdiff(1:n,ec);
t=cell(length(r),1);
nt=[];
nn=[];
nc=[];
for i=r
  tpos=regexp(colnames{i},'\(text\)');
  npos=regexp(colnames{i},'\(numeric\)');
  if ~isempty(tpos)
    colnames{i}=colnames{i}(1:(tpos-1));
    t{i}='text';
    nt=[nt i];
  elseif ~isempty(npos)
    colnames{i}=colnames{i}(1:(npos-1));
    t{i}='numeric';
    nn=[nn i];
  else
    t{i}='categoric';
    nc=[nc i];
  end
end

for j=1:length(nt)
  si.supacc{j}=colnames{nt(j)};
  si.supdesc{j}=colnames{nt(j)};
end

first=get_cell_elem(f,1);
first=cat(1,first{:});
er=find(cellfun('isempty',first));
ur=setdiff(1:length(f),er);
f=f(ur);

for j=1:length(nt)
  si.supdat{j}=get_cell_elem(f,nt(j))';
end

for i=length(ur)
  if ~isempty(f{i}{1})
      si.supdat{j,cs}=f{i}{nt(j)};
    end
    
    cs=cs+1;
  end
end

