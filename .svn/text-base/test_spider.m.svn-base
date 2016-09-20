clear;
tts=[2     8    10    13    22];
tts=10;
for i=1:length(tts)
  tt=tts(i);
  load ~/projects/miRNA/new/GCM/miR_GCM_237_th32_log2_D.mat
  Dtrain=D;
  tt_supid=strmatch('TT',Dtrain.supacc);
  out=double(Dtrain.supdat(tt_supid,:)==tt)';
  Dtrain.spider.dat=data('train',Dtrain.dat',2*(out-0.5));
  
  load ~/projects/miRNA/new/PD/miR_PD_19_th32_log2_D.mat
  Dtest=D;
  tt_supid=strmatch('TT',Dtest.supacc);
  out=double(Dtest.supdat(tt_supid,:)==tt)';
  Dtest.spider.dat=data('test',Dtest.dat',2*(out-0.5));
  
  [r,a]=train(knn({distance('euclid'),'k=3'}),Dtrain.spider.dat);
  loss(r)
  res=test(a,Dtest.spider.dat);
  tt_supid=strmatch('TT',D.supacc);
  cls0=find(D.supdat(tt_supid==1));
  
  
  X=Dtrain.dat';
  nX=size(X,1);
  XO=double(Dtrain.supdat(tt_supid,:)==tt)';
  Y=Dtest.dat';
  YO=double(Dtest.supdat(tt_supid,:)==tt)';
  nY=size(Y,1);
  
  c=knn_classify(Y,X,[XO'; 1-XO'],ones(2,1),3);
  % [c,s]=knn_classify(X,X,[XO'; 1-XO'],ones(2,1),3);
  
  disp(i);
  find(c==1)
  find(get_x(res)>0)'
  loss(res)
  crosstab(get_x(res)>0,get_y(res)>0)
end


XY=[X; Y];
XYO=[XO; YO];
d=dna_dist(XY);
d=d(1:nX,nX+(1:nY));

[ds,idx]=sort(d);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test cv
clear

load ~/projects/miRNA/new/GCM/miR_GCM_237_th32_log2_D.mat
Dtrain=D;
tt_supid=strmatch('TT',Dtrain.supacc);
tt=8;
Dtrain=reorder_D_cols(Dtrain,intersect(setdiff(1:size(Dtrain.dat,2),...
                                               find(cat(1,Dtrain.sup(:).PD))),...
                                       find((cat(1,Dtrain.sup(:).MAL)==2) ...
                                            & (cat(1,Dtrain.sup(:).CLT)==1))));


a1=param(knn({distance('euclid')}),'k',1:2:11);
a2=gridsel(a1);
[r,a]=train( cv(knn({distance('euclid'),'k=3'}),'folds=5'), ...
             data('train',Dtrain.dat',2*((Dtrain.supdat(tt_supid,:)==tt)'-0.5)));
get_mean(r)
