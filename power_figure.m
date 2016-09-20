function [P,C,K,P1,C1,K1]=power_figure(freq,fig_params)

clf;
if ~exist('fig_params','var') || ~isfield(fig_params,'type')
  fig_params.type=2;
end


default.r0=1.2e-6*1500*3;
default.Ng=13023;
default.rnk=6;
default.ND=50;
default.Dstep=1;
default.Dplot=5;
default.NV=200;
default.Vstep=2;
default.Vplot=5;
default.L=1500;
default.q=0.1;
default.norm_frac1=0.125;
default.norm_frac2=1;
default.cost_per_amp=2;
default.fnr=0.35;
default.pointD=11;
default.pointV=24;
default.ext='';

fig_params=add_struct(default,fig_params);
ND=fig_params.ND;
Dstep=fig_params.Dstep;
Dplot=fig_params.Dplot;
NV=fig_params.NV;
Vstep=fig_params.Vstep;
Vplot=fig_params.Vplot;
rnk=fig_params.rnk;
Ng=fig_params.Ng;
ext=fig_params.ext;
L=fig_params.L;
q=fig_params.q;
norm_frac1=fig_params.norm_frac1;
norm_frac2=fig_params.norm_frac2;
cost_per_amp=fig_params.cost_per_amp;
fnr=fig_params.fnr;
pointD=fig_params.pointD;
pointV=fig_params.pointV;

disp(fig_params);

r0=fig_params.r0;
r1=r0+(freq*(1-fnr));

P=zeros(ND/Dstep,NV/Vstep);
P1=P;
C=zeros(ND/Dstep,NV/Vstep);
C1=C;
K=zeros(ND/Dstep,NV/Vstep);
K1=K;
for i=1:ND/Dstep
  disp([ num2str(i) '/' num2str(ND/Dstep)]);
  for j=1:NV/Vstep
    n1=i*Dstep;
    n2=j*Vstep;
    [P1(i,j),K1(i,j)]=single_screen_power_smooth(r0,r1,n1+n2,Ng,rnk,q,L);
    [P(i,j),K(i,j)]=screen_power_smooth(r0,r1,n1,n2,Ng,rnk,q,L,K1(i,j),0.1);
    C1(i,j)=screen_cost(r0,n1+n2,0,Ng,norm_frac2,0,15,cost_per_amp); % we use norm_frac2 for single phase
    C(i,j)=screen_cost(r0,n1,n2,Ng,norm_frac1,norm_frac2,15,cost_per_amp);
  end
end
switch fig_params.type
 case {1,2},
  imagesc(P);
 case {3,4},
  surf(P);
  axis([1 NV/Vstep 1 ND/Dstep 0 1]);
  shading interp
  hidden on
end
caxis([0 1]);
title(['Power to detect ' num2str(freq)]);
ylabel('N Discovery');
xlabel('N Validation');
set(gca,'YTick',1:Dplot:ND/Dstep,'YTickLabel',cellstr(num2str((1:Dplot:ND/Dstep)'*Dstep)));
set(gca,'XTick',1:Vplot:NV/Vstep,'XTickLabel',cellstr(num2str((1:Vplot:NV/Vstep)'*Vstep)));
set(gca,'YDir','normal');
cb=colorbar('SouthOutside');
pos_cb=get(cb,'Position');
pos1=get(gca,'Position');
set(cb,'Position',[ pos_cb(1) pos_cb(2)-0.03 pos_cb(3) 0.03]);
set(gca,'Position',pos1);

hold on;
set(gcf,'renderer','painters');


% surf(1:NV/Vstep,1:ND/Dstep,P,C);

c=contourc(C,15);
st=1;
for k=1:15
  disp(c(1,st))
  disp([ st c(2,st)]);
  po=zeros(1,size(c,2));   
  po1=po;
  ko1=po;
  for kk=(st+1):(st+c(2,st))
    [po1(kk),ko1(kk)]=single_screen_power_smooth(r0,r1,round(c(2,kk)*Dstep)+round(c(1,kk)*Vstep),Ng,rnk,q,L);
    po(kk)=screen_power_smooth(r0,r1,round(c(2,kk)*Dstep),round(c(1,kk)*Vstep),Ng,rnk,q,L,ko1(kk));
  end
  [mpo,mpi]=max(po((st+1):(st+c(2,st))));
  switch fig_params.type
   case 2
    lh=line(c(1,(st+1):(st+c(2,st))),c(2,(st+1):(st+c(2,st))));
    set(lh,'Color',[0 0 0],'LineWidth',1);
    if isfield(fig_params,'max_power') && fig_params.max_power==1
      plot(c(1,st+mpi),c(2,st+mpi),'wx');
    end
    ax=axis;
    if (0)
      if c(1,st)<1e7
        text(ax(2)+(ax(2)-ax(1))*0.001,c(2,st+c(2,st)),[sprintf('$%2.1f',c(1,st)/1000) 'K'],'FontSize',8);
      else
        text(ax(2)+(ax(2)-ax(1))*0.001,c(2,st+c(2,st)),[sprintf('$%2.1f',c(1,st)/1e6) 'M'],'FontSize',8);
      end     
    end
   case 4 
    lh=line(c(1,(st+1):(st+c(2,st))),c(2,(st+1):(st+c(2,st))),po((st+1):(st+c(2,st)))+0.001);
    set(lh,'Color',[0 0 0],'LineWidth',1);
    plot3(c(1,st+mpi),c(2,st+mpi),po(st+mpi)+0.001,'wx');
  end
  st=st+c(2,st)+1;
end

if (1) % show V&V point
  ax=axis;
  if length(ax)==4
    ph=plot(pointV/Vstep,pointD/Dstep,'w.');
    if isfield(fig_params,'point_power') && fig_params.point_power==1
      th=text((pointV+1)/Vstep,(pointD+1)/Dstep,num2str(P(pointD/Dstep,pointV/Vstep)),'Color',[1 1 1]);
    end
  else
    ph=plot3(24/Vstep,11/Dstep,screen_power_smooth(r0,r1,pointD,pointV,Ng,rnk,q,L)+0.001,'w.');
    if isfield(fig_params,'point_power') && fig_params.point_power==1
      th=text(24/Vstep,11/Dstep,screen_power_smooth(r0,r1,pointD,pointV,Ng,rnk,q,L)+0.1,num2str(screen_power_smooth(r0,r1, ...
                                                        pointD,pointV,Ng,rnk,q,L)),'Color',[1 1 1]);
    end
  end
end

print_D(['power_' num2str(freq) ext],{{'fig'},{'pdf'},{'png','-r180'}});




