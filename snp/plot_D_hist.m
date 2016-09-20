function plot_D_hist(CL21,fname,x_axis)
CL21=shrink_D(CL21);

perpage_x=4;
perpage_y=6;
perpage=perpage_x*perpage_y;
close all
gr=make_subplotgrid(4*ones(1,perpage_x),3*ones(1,perpage_y),...
                    0.5*ones(1,perpage_x+1),1.2*ones(1,perpage_y+1),0.3,0.3);
if exist([fname '.ps'],'file');
  delete([fname '.ps']);
end
if exist([fname '.pdf'],'file');
  delete([fname '.pdf']);
end

for i=1:perpage:size(CL21.dat,2)
  figure(1); clf;
  for ii=i:min(i+perpage-1,size(CL21.dat,2))
    subplotgrid(gr,floor((ii-i)/perpage_x)+1,mod(ii-i,perpage_x)+1);
    hist(CL21.dat(:,ii),50);
    if exist('x_axis','var')
      ax=axis;
      axis([ x_axis ax(3:4)]);
    end
    title(['Sample ' num2str(ii)]);
  end
  rend=get(gcf,'renderer');
  set(gcf,'renderer','none');
  print('-f1','-dpsc2','-append',[ fname '.ps']);
  set(gcf,'renderer',rend);  
%  print_D([ figdir 'smoothed_' num2str(i) '_' num2str(min(i+perpage-1,size(CL21.dat,2)))],...
%          {{'pdf'}});
end
ps2pdf([fname '.ps'],[fname '.pdf'],1);
