function [ph,th]=plot_pca(dat,markset,labelset,k)
if nargin<4
  k=2;
end

[prj,m,D,V]=pca(dat,k);  
if ( prj==-1) % too few eigenvalues
  return;
end
if (size(dat,1)>1000) 
  perm=randperm(size(dat,1));
  slct_set=perm(1:1000);
end
hold off
if nargin>2 & ~isempty(labelset)
  if k==2
    ph=plot(prj(:,1),prj(:,2),'.','MarkerSize',1);
  else
    ph=plot3(prj(:,1),prj(:,2),prj(:,3),'.','MarkerSize',1);
  end
  hold on
  for i=1:size(labelset,1)
    if k==2
      th(i)=text(prj(i,1),prj(i,2),labelset{i});
    else
      th(i)=text(prj(i,1),prj(i,2),prj(i,3),labelset{i});
    end
  end
else
  if k==2
    ph=plot(prj(:,1),prj(:,2),'x')
  else
    ph=plot3(prj(:,1),prj(:,2),prj(:,3),'x')  
  end
  hold on
  
  if nargin>1 & ~isempty(markset)
    if ~iscell(markset)
      markerset={markset};
    end
    marks={'ro','bs','kv','c^'};
    if k==2
      for i=1:length(markset)
        mh{i}=plot(prj(markset{i},1),prj(markset{i},2),marks{mod(i-1, ...
                                                          length(marks))+1});
      end
    else
      for i=1:length(markset)
        mh{i}=plot3(prj(markset{i},1),prj(markset{i},2),prj(markset{i},3),marks{mod(i-1, ...
                                                          length(marks))+1});
      end
    end
  end
end
title(['First ' num2str(k) ' PCA components'])
axis equal
% enlarge_axis(0.1);
