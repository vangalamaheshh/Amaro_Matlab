function res = draw_nmf_radial_plot(res,P)
% based on mut/analysis/20110909_pancan/run.m

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'radial_plot_fontsize',3);
P=impose_default_value(P,'draw_nmf_radial_plot_method',2);

demand_field(res,{'pat','patfrac','ttype'});
demand_field(res.pat,{'ntot','ttype','ttype_idx'});
if ~isfield(res.ttype,'color')
  res.ttype.color = distinct_colors(slength(res.ttype));
end

nf = size(res.patfrac,2);

res.rho = -0.6+max(0,min(4,log10(res.pat.ntot)));
frac = res.patfrac;

% map to circumference

method = P.draw_nmf_radial_plot_method;
if method==1

  z = 6*ones(size(frac));
  z(frac>0.1) = 5;
  z(frac>0.2) = 4;
  z(frac>0.3) = 3;
  z(frac>0.5) = 2;
  z(frac>0.8) = 1;
  [tmp ord] = sort(z,2);
  res.theta = ((ord(:,1)-1)/nf)*2*pi + ...
          ((ord(:,2)-1)/nf)*(0.9*2*pi/nf) + ...
          ((ord(:,3)-1)/nf)*(0.9*2*pi/nf/nf) + ...
          ((ord(:,4)-1)/nf)*(0.9*2*pi/nf/nf/nf) + ...
          ((ord(:,5)-1)/nf)*(0.9*2*pi/nf/nf/nf/nf);
  
elseif method==2

  z = frac;
  z = round(z*20)/20;
  z(z<0.1)=0;
  [tmp ord] = sort(z,2,'descend');
  [u ui uj] = unique(ord,'rows');
  res.theta = ((uj-1)/length(u))*2*pi;

elseif method==3

  
  


end

[res.x res.y] = polar_to_cartesian(res.rho,res.theta);

a.old = { 'AML' 'Bladder' 'Breast' 'CLL' 'CRC' 'Carcinoid' 'Cervical' 'DLBCL' 'Esoph' 'Ewing',...
          'GBM' 'HeadNeck' 'KidneyRC' 'KidneyRP' 'LGG' 'LUAD' 'LUSC' 'MM' 'Medullo' 'Melanoma',...
          'NB' 'OV' 'Pancreas' 'Prostate' 'Rhabdoid' 'Stomach' 'Thyroid' 'Uterine'};
a.new = { 'AML' 'Blad' 'Bre' 'CLL' 'CRC' 'Car' 'Cerv' 'DLB' 'Eso' 'Ew',...
          'GBM' 'HN' 'KRC' 'KRP' 'LGG' 'LUAD' 'LUSC' 'MM' 'Med' 'Mel',...
          'NB' 'Ov' 'Panc' 'Pro' 'Rh' 'Sto' 'Thy' 'Ute'};
labels = apply_aliases(res.pat.ttype,a);

clf;hold on
for i=1:slength(res.ttype)

  
  idx=find(res.pat.ttype_idx==i);
  scatter(res.x(idx),res.y(idx),5,res.ttype.color(i,:),'filled');
  %    scatter(x(idx),y(idx),20,res.ttype.color(i,:),'filled');
%  idx = find(res.pat.ttype_idx==i);
% text(res.x(idx),res.y(idx),labels(idx),'color',res.ttype.color(i,:),'fontsize',P.radial_plot_fontsize);

end
hold off,set(gca,'visible','off');
legend(res.ttype.name,'location','eastoutside'), legend('boxoff');
set(gcf,'color',[1 1 1],'position',[292 25 1130 793]);

xlim([min(res.x) max(res.x)]);ylim([min(res.y) max(res.y)]);





