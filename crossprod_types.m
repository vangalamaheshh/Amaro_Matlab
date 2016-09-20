function [D,cross_supid]=crossprod_types(D,t1,t2)

s=crosstab_sets(D.supdat(t1,:),D.supdat(t2,:));

if is_multiple_types(D,t1)
  [typeacc1,typedesc1,D1,r1]=decollapse_supdat(D,t1);
  acc1=D1.supacc(r1,:);
  desc1=D1.supdesc(r1,:);
else
  typeacc1=deblank(D.supacc(t1,:));
  typedesc1=deblank(D.supacc(t1,:));
  acc1=strvcat(typeacc1,['n' typeacc1]);
  desc1=strvcat(typedesc1,['Not ' typedesc1]);
end

if is_multiple_types(D,t2)
  [typeacc2,typedesc2,D2,r2]=decollapse_supdat(D,t2);
  acc2=D2.supacc(r2,:);
  desc2=D2.supdesc(r2,:);
else
  typeacc2=deblank(D.supacc(t2,:));
  typedesc2=deblank(D.supacc(t2,:));
  acc2=strvcat(typeacc2,['n' typeacc2]);
  desc2=strvcat(typedesc2,['Not ' typedesc2]);
end

v=NaN*ones(1,size(D.dat,2));
[fi,fj]=find(cell2mat(cellfun_any('~isempty',s)));
combacc=[deblank(typeacc1) 'x' deblank(typeacc2) ': '];
combdesc=[deblank(typedesc1) ' x ' deblank(typedesc2) ': '];
for i=1:length(fi)
  combacc=[combacc [deblank(acc1(fi(i),:)) '-' ...
                    deblank(acc2(fj(i),:))]  ];
  combdesc=[combdesc [deblank(desc1(fi(i),:)) '-' ...
                    deblank(desc2(fj(i),:))]  ];
  if (i<length(fi))
    combacc=[ combacc '/'];
    combdesc=[ combdesc '/'];
  end    
  v(s{fi(i),fj(i)})=i;
end

D.supacc=strvcat(D.supacc,combacc);
D.supdesc=strvcat(D.supdesc,combdesc);
D.supdat=[ D.supdat; v ];
cross_supid=size(D.supdat,1);
% decollapse the types and the crossprod


