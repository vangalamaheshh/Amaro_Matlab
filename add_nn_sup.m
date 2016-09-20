function [Dtest,range]=add_nn_sup(Dtrain,Dtest,supid,nnbr,dist_type,fet_sel_supid,nfet,selection_params,xv_flag)

knn=zeros(nnbr,size(Dtest.dat,2));
if ~isempty(selection_params)
   si=gg_sel_fet(Dtrain,fet_sel_supid,nfet,selection_params);
%    s(si)
%    Dtrain.gacc(si,:)
   disp('adding nn sup using');
   disp(strvcat(Dtrain.gacc(si)));
else
    si=1:size(Dtrain.dat,1);
end

for i=1:size(Dtest.dat,2)
  if xv_flag
    subset=setdiff(1:size(Dtrain.dat,2),i);
    if ~isempty(selection_params)
      si=gg_sel_fet(reorder_D_cols(Dtrain,subset), ...
                    fet_sel_supid,nfet,selection_params);
    end
  else
    subset=1:size(Dtrain.dat,2);
  end
  dat=Dtrain.dat(si,subset);
  
  d=dist(Dtest.dat(si,i)',dat',dist_type);
  [ds,dord]=sort(d);
  knn(:,i)=subset(dord(1:nnbr));
end 

range=[];
for i=1:nnbr
    [Dtest,supidx]=add_D_sup(Dtest,['NBR' num2str(i)],['Neighbor #' ...
                        num2str(i)],Dtrain.supdat(supid,knn(i,:)),'col');
    range=[range supidx];
end
