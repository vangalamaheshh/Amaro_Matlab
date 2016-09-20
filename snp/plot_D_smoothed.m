function plot_D_smoothed(CL21,cyto)
CL21=shrink_D(CL21);

perpage_x=4;
perpage_y=6;
perpage=perpage_x*perpage_y;
close all
gr=make_subplotgrid(4*ones(1,perpage_x),3*ones(1,perpage_y),...
                    0.5*ones(1,perpage_x+1),1.2*ones(1,perpage_y+1),0.3,0.3);
for i=1:perpage:size(CL21.dat,2)
  figure(1); clf;
  for ii=i:min(i+perpage-1,size(CL21.dat,2))
    subplotgrid(gr,floor((ii-i)/perpage_x)+1,mod(ii-i,perpage_x)+1);
    plot_smoothed(reorder_D_cols(CL21,ii),1,cyto);
    ii
  end
  rend=get(gcf,'renderer');
  set(gcf,'renderer','none');
  print('-f1','-dpsc2','-append','CL21_smoothed.ps');
  set(gcf,'renderer',rend);  
%  print_D([ figdir 'smoothed_' num2str(i) '_' num2str(min(i+perpage-1,size(CL21.dat,2)))],...
%          {{'pdf'}});
end
ps2pdf('CL21_smoothed.ps','CL21_smoothed.pdf',1);
