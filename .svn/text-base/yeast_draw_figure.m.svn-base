function yeast_draw_figure(Q,lab,true_labs,labelled_pts,unlabelled_pts,J,q,has_labels,figname)

figure(1); clf; set(gcf,'renderer','zbuffer');

h1=subplot(1,10,1:9);
[sind1,sind1i]=sort(lab(unlabelled_pts));
svec2=[labelled_pts' unlabelled_pts(sind1i)];
Q2=Q(svec2,svec2);
spy(Q2); 
axis square;
hold on;
Q3=sparse(size(Q2,1),size(Q2,2));
Q3(1:length(labelled_pts),1:length(labelled_pts))=Q2(1:length(labelled_pts),1:length(labelled_pts));
spy(Q3,'r');

Q4=sparse(size(Q2,1),size(Q2,2));
for i=1:max(lab)
  pred_as_i=find(lab(svec2)==i);  
  Q4(pred_as_i,pred_as_i)=Q2(pred_as_i,pred_as_i);
end
spy(Q4,'k');
xlabel(['Black - within type, Red - between types in train, Blue - ' ...
        'between types in test' ]);

h2=subplot(1,10,10);
L=[lab true_labs lab==true_labs];
%sum(L)
imagesc(L(svec2,:));
caxis([0 max(max(L))]);
cm=colormap;
cm(1,:)=[1 1 1];
colormap(cm);
text(1,size(L,1)*1.01,'Assignment','Rotation',270);
text(2,size(L,1)*1.01,'True','Rotation',270);
text(3,size(L,1)*1.01,'Correct','Rotation',270);
set(gca,'XTick',[]);

axes(h1);

e=calc_configuration_energy([1:q lab'],J,q,has_labels);
title([figname ' energy=' num2str(e) ' number of correct=' num2str(length(find(lab==true_labs))-has_labels+1)]);
