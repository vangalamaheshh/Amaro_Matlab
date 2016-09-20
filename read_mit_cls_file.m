function D=read_mit_cls_file(D,fname,all_in_one_sup)

if ~exist('all_in_one_sup','var')
  all_in_one_sup=0;
end

f=read_dlm_file(fname,' ');

n=str2num(f{1}{1});
nt=str2num(f{1}{2});
if f{2}{1}(1)=='#' 
    has_names=1;
    for i=1:nt
        nms{i}=f{2}{i+1};
    end
else
    has_names=0;
end 

for i=1:n
  cls(i)=str2num(f{2+has_names}{i});
end

if has_names
  if all_in_one_sup
    D.supacc='PHEN: ';
    D.supdesc='Phenotype: ';
    for i=1:(length(nms)-1)
      D.supacc=[ D.supacc num2str(i) '-' nms{i} '/' ];
      D.supdesc=[ D.supdesc num2str(i) '-' nms{i} '/' ];
    end
    D.supacc=[ D.supacc num2str(length(nms)) '-' nms{end} ];
    D.supdesc=[ D.supdesc num2str(length(nms)) '-' nms{end} ];  
  else
    D.supdesc=strvcat(nms);
    D.supacc=strvcat(nms);
  end
else
    D.supacc=[repmat('CLS',n,1) num2str(cls')];
    D.supdesc=[repmat('Class ',n,1) num2str(cls')];
end

if all_in_one_sup
  D.supdat=cls;
  if min(cls)==0
    D.supdat=D.supdat+1;
  end
else
  D.supdat=zeros(nt,n);
  for i=1:n
    D.supdat(cls(i)+1,i)=1;
  end
end

