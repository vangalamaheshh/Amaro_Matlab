function plot_std_vs_mean(D,use_median,supid)

if ~exist('supid','var');
  supid=0;
end

if ~exist('use_median','var');
  use_median=0;
end

if supid>0
  u=unique(D.supdat(supid,:));
  for i=1:length(u)
    plot_std_vs_mean(reorder_D_cols(D,find(D.supdat(supid,:)==u(i))),use_median,0);
    hold on;
  end
else
  if use_median
    s=mad(D.dat,1,2);
    m=median(D.dat,2);
    plot(m,s,'.');
    xlabel('Median');
    ylabel('MAD');
  else
    s=std(D.dat,0,2);
    m=mean(D.dat,2);    
    plot(m,s,'.');
    xlabel('Std.');
    ylabel('Mean');
  end
end
