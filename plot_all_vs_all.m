function plot_all_vs_all(x,labels)

n=size(x,1);
n
for i=1:(n-1)
  for j=(i+1):n
    subplot(n,n,j+(i-1)*n);
    plot(x(j,:),x(i,:),'.');
    if exist('labels','var')
      xlabel(labels{j});
      ylabel(labels{i});
    else
      xlabel(num2str(j));
      ylabel(num2str(i));
    end
  end
end
