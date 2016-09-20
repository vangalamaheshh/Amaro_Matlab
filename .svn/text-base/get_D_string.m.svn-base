function str=get_D_string(D,idx,st)
% assumes for cols

Dx=D;
n=regexp(st,'\$\{(\w*)\}*','tokens');
for i=1:length(n)
  sacc=n{i}{1};
  si=strmatch([sacc ':'],D.supacc);
  [typeacc,typedesc,Dx,range,non_empty]=decollapse_supdat(D,si);
  [mv,mi]=max(Dx.supdat(range,idx),[],1);
  for j=1:length(idx)
    if mv==0
      s{i,j}='EMPTY';
    else
      s{i,j}=deblank(Dx.supdesc(range(mi(j)),:));
    end
  end
  tok{i}=['\$\{' sacc '\}'];
end

for i=1:length(idx)
  str{i}=regexprep(st,tok,s(:,i));
end


