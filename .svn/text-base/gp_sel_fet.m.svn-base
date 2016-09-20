function fetid=gp_sel_fet(D,supid,nfet,selection_params)


if ischar(selection_params)
  selection_params.method='gp';
  selection_params.score=selection_params;
end

switch selection_params.method
 case 'gp'
  cls=unique_keepord(D.supdat(supid,:));
  ncls=length(cls);
  
  fetpercls=get_parts(1:nfet,ncls);
  fetid=zeros(nfet,1);
  for i=1:ncls
    [p,s]=differential_analysis(D,find(D.supdat(supid,:)~=cls(i)),find(D.supdat(supid,:)==cls(i)),'gcsnr');
    [ss,si]=sort(s);
    si=flipud(si);
    fetid(fetpercls{i})=si(1:length(fetpercls{i}));
  end
 case '2class'
  % assume only 2 classes
  [p,s]=differential_analysis(D,find(D.supdat(supid,:)==0),find(D.supdat(supid,:)==1),'gcsnr');
  [ss,si]=sort(abs(s));
  si=flipud(si);
  fetid=si(1:nfet);
end
