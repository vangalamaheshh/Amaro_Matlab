function [ap2,p2,p1_gte,p1_lte,sets,s,n,acc,select_compounds1,select_compounds2]=...
    generate_data_for_figure(D,cls0,cls1,fname,do_perm)
if (0)
    [ap2,s]=differential_analysis(D,cls0,cls1,struct('method','ttest_minvar','minvar',0.4^2),0);
    p2=[];
    p1_gte=[];
    p1_lte=[];
    sets=[];
    n=[];
    acc=[];
    select_compounds1=[];
    select_compounds2=[];
    return
end

if ~exist('do_perm','var')
  do_perm=1;
end

if ~do_perm  
  [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr,topidx]=...
      marker_selection(D.dat,cls0,cls1,...
                       struct('method','ttest_minvar','minvar',0.4^2,...
                              'nparts_perm',1,'nparts_fix',0),-1);
    ap2=Pf;
    p1_gte=P1sgte;
    p1_lte=P1slte;
    p2=P2s(:,1);
    s=S;
    acc=[];
    sets=[];
    n=[];
    
    select_compounds1=find(p2<0.0002*2 & s<0 );
    [dum,ord]=sort(s(select_compounds1));
    select_compounds1=select_compounds1(ord);
    
    select_compounds2=find(p2<0.0002*2 & s>0);
    [dum,ord]=sort(s(select_compounds2));
    select_compounds2=select_compounds2(flipud(ord));
    keyboard
    
    save([fname '.f.mat'],'p1_gte','p1_lte','p2','ap2','sets','s','n','acc','select_compounds1','select_compounds2');
  
else

  sets{1}=1:size(D.dat,1);
  p2=zeros(size(D.dat,2));
  n=[500 5000 50000 200000];
  for i=1:2
    %1 all x 500
    if i>1
      sets{i}=find(p2<10/n(i-1));
    end 
    D1=reorder_D_rows(D,sets{i});
    disp(['running ' num2str(n(i)) ' permutations on ' num2str(length(sets{i})) ' features']);
    [P1sgte,P1slte,P2s,Pf,rs,Pf2,S,gp,fwer,fpr,topidx]=...
        marker_selection(D1.dat,cls0,cls1,...
                         struct('method','ttest_minvar','minvar',0.4^2,...
                                'nparts_perm',10,'nparts_fix',10),n(i));
    ap2(sets{i})=Pf(:,1);
    p1_gte(sets{i})=P1sgte(:,1);
    p1_lte(sets{i})=P1slte(:,1);
    p2(sets{i})=P2s(:,1);
    s(sets{i})=S;
    acc(sets{i})=i;
  end

  select_compounds1=find(p2<0.0002*2 & p1_lte<0.5 & ~ ...
                         isnan(p2));
  %[dum,ord]=sort(mdiff(topidx(select_compounds1)));
  [dum,ord]=sort(s(select_compounds1));
  select_compounds1=select_compounds1(ord);

  select_compounds2=find(p2<0.0002*2 & p1_lte>0.5 & ~isnan(p2) );
  %[dum,ord]=sort(mdiff(select_compounds2));
  [dum,ord]=sort(s(select_compounds2));
  select_compounds2=select_compounds2(flipud(ord));

  save([fname '.f.mat'],'p1_gte','p1_lte','p2','ap2','sets','s','n','acc','select_compounds1','select_compounds2');
end


