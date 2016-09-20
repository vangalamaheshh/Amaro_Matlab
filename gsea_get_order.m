function [ord,pv,st]=gsea_get_order(dat,phen,test)

if ischar(test)
  test.method=test;
end

D.dat=dat;

[P,S]=differential_analysis(D,find(~(phen==1)),find(phen==1),test,0);
[tmp,ord]=sort(P);

return

%%% PREV VERSION

switch(lower(test))
 case 'ttest'
  [pv,st]=ttest2_many(dat,find(phen==1),find(~(phen==1)));
  [tmp,ord]=sort(st);
  ord=flipud(ord);
 case 'tnom'
  [pv,nerr]=ttest2_many(dat,find(phen==1),find(~(phen==1)));
  [tmp,ord]=sort(nerr);
 case 'ttest2'
  [pv,st]=ttest2_many(dat,find(phen==1),find(~(phen==1)));
  [tmp,ord]=sort(abs(st));  
  ord=flipud(ord);
 case 'permttest2'
  [Pf,rs,Pf2,S,gp,fwer,fpr]=marker_selection(dat,find(phen==1),find(~(phen==1)),'ttest',1000);
  [tmp,ord]=sortrows([Pf S]);
  ord=flipud(ord);  
end
