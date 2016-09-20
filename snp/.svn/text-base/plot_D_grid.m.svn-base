function plot_D_grid(fname,CL21,fnc,perpage_x,perpage_y,varargin)
% CL21=shrink_D(CL21);

if ~exist('perpage_x','var') || isempty(perpage_x) || perpage_x==0
  perpage_x=4;
end

if ~exist('perpage_y','var') || isempty(perpage_y) || perpage_y==0
  perpage_y=6;
end

perpage=perpage_x*perpage_y;
close all
gr=make_subplotgrid(4*ones(1,perpage_x),3*ones(1,perpage_y),...
                    0.5*ones(1,perpage_x+1),1.2*ones(1,perpage_y+1),0.3,0.3);
for i=1:perpage:size(CL21.dat,2)
  figure(1); clf;
  for ii=i:min(i+perpage-1,size(CL21.dat,2))
    ii
    subplotgrid(gr,floor((ii-i)/perpage_x)+1,mod(ii-i,perpage_x)+1);
    if exist('varargin','var')
      fnc(reorder_D_cols(CL21,ii),varargin{:});
    else
      fnc(reorder_D_cols(CL21,ii));
    end
  end
  rend=get(gcf,'renderer');
  set(gcf,'renderer','none');
  print('-f1','-dpsc2','-append',[ fname '.ps']);
  set(gcf,'renderer',rend);  
%  print_D([ figdir 'smoothed_' num2str(i) '_' num2str(min(i+perpage-1,size(CL21.dat,2)))],...
%          {{'pdf'}});
end
ps2pdf([ fname '.ps'],[ fname '.pdf'],1);
