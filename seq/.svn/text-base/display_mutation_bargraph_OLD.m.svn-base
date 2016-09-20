function display_mutation_bargraph(D,P)

if ~exist('P','var'), P=[]; end
P=impose_default_value(P,'by','rate');
P=impose_default_value(P,'omit_silent',false);
P=impose_default_value(P,'ymax',[]);
P=impose_default_value(P,'legend','best');

names = D.names;
nnon_tot = D.nnon_tot;
nsil_tot = D.nsil_tot;
rate_sil = D.rate_sil;
rate_non = D.rate_non;
N_tot = D.N_tot;

omit_legend = P.omit_silent;

if strcmp(P.by,'count')
  if P.omit_silent
    bar(nnon_tot);
    y = nnon_tot;
  else
    bar([nnon_tot nsil_tot],'stacked')
    y = nnon_tot+nsil_tot;
  end
  title('mutation counts','fontsize',20);
  ylabel('number of mutations','fontsize',15);

elseif strcmp(P.by,'rate')
  if P.omit_silent
    bar(1e6*rate_non);
    y = 1e6*rate_non;
  else
    bar(1e6*[rate_non rate_sil],'stacked')
    y = 1e6*(rate_non+rate_sil);
  end
  title('mutation rates','fontsize',20);
  ylabel('mutations per total Mb','fontsize',15);

elseif strcmp(P.by,'coverage')
  bar(N_tot/1e6);
  y = N_tot/1e6;
  title('sequencing coverage','fontsize',20);
  ylabel('total Mb sequenced','fontsize',15);
  omit_legend = true;
else
  error('unknown P.by');
end

if ~omit_legend
  legend('nonsilent','silent','location',P.legend);
end

xlim([0.2 length(names)+0.8]);
set(gca,'xtick',1:length(names),'xticklabel',names,'fontsize',10);
set(gca,'TickLength',[0 0]);
xticklabel_rotate(0.33+(1:length(names)),90,names,'fontsize',10);

if ~isempty(P.ymax)
  ylim([0 P.ymax]);
  for i=1:length(names)
    if y(i)>P.ymax
      text(i+0.2,P.ymax*0.95,num2str(round(y(i))),'horizontalalignment','center','rotation',90,...
            'verticalalignment','middle','color',[1 1 1]);
    end
  end
end

set(gca,'position',[0.05 0.1 0.9 0.8])
