function [typeacc,typedesc,D,range,non_empty]=decollapse_supdat(D,si,force_vals,add_sup_field)

posacc=find(D.supacc(si,:)==':');
if isempty(posacc)
  typeacc=D.supacc(si,:);
  typedesc=D.supdesc(si,:);
  range=si;
  non_empty=1;
  return;
end
% assert(~isempty(posacc));
typeacc=D.supacc(si,1:(posacc(1)-1));

posdesc=find(D.supdesc(si,:)==':'); 
assert(~isempty(posdesc));
typedesc=D.supdesc(si,1:(posdesc(1)-1));

sacc=posacc(1)+2;
pacc=find([ deblank(D.supacc(si,:)) '/']=='/');
assert(~isempty(pacc));
sdesc=posdesc(1)+2;
pdesc=find([deblank(D.supdesc(si,:)) '/']=='/');
assert(~isempty(pdesc));
assert(length(pacc)==length(pdesc));

if exist('force_vals','var') && length(force_vals)==1 && force_vals==-1
  force_vals=1:(length(pacc)+1);
end

range=[];
non_empty=[];
for i=1:length(pacc)
  curacc=D.supacc(si,sacc:(pacc(i)-1));
  dashpos=find(curacc=='-');
  if ~isempty(dashpos)
    curacc=curacc((dashpos(1)+1):end);
  end
  curdesc=D.supdesc(si,sdesc:(pdesc(i)-1));
  dashpos=find(curdesc=='-');
  if ~isempty(dashpos)
    curdesc=curdesc((dashpos(1)+1):end);
  end  
%  disp([ num2str(i) ' ' curacc ' ' curdesc '.']);
  sacc=pacc(i)+1;
  sdesc=pdesc(i)+1;
  v=D.supdat(si,:)==i;
  if ( (~exist('force_vals','var') && (sum(v)>0)) | ...
       (exist('force_vals','var') && ismember(i,force_vals)) )
    D.supdat=[ D.supdat; D.supdat(si,:)==i ];
    non_empty=[ non_empty; i];
    range=[ range; size(D.supdat,1)];
    D.supacc=strvcat(D.supacc,curacc);
    D.supdesc=strvcat(D.supdesc,curdesc);
  else
    verbose(['No cases from type ' curacc ' [' curdesc ']'],2);
  end
end


sdm=mat2cell(D.supdat,ones(size(D.supdat,1),1),ones(size(D.supdat, ...
                                                  2),1));
for i=1:size(D.supacc)
  pos=find(D.supacc(i,:)==':');
  if ~isempty(pos)
    spacc{i}=deblank(D.supacc(i,1:(pos(1)-1)));
  else
    spacc{i}=deblank(D.supacc(i,:));    
  end
  tmp=spacc{i};
  tmp(tmp=='-')='_';
  tmp=regexprep(tmp,'+','_PLUS_');
  if ismember(tmp(1),'0':'9')
      tmp=['X' tmp];
  end   
  spacc{i}=tmp;
end

if exist('add_sup_field','var') && add_sup_field
  D.sup=cell2struct(sdm,regexprep(regexprep(spacc,'[^0-9a-zA-Z]','_'),'^([^a-zA-Z])','X$1'),1);
end







