function bsub_log_hist
x = load_lines('/xchip/tcga_scratch/lawrence/log/bsub.log');
t = parse(x,'^trying to mount (\d*)$',{'num'});
j = t.num;
[u ui uj] = unique(j);
for i=1:length(u), uc(i,1) = sum(uj==i); end
hist(uc,min(uc):max(uc));
title('attempts to mount /xchip/cga1','fontsize',15);
xlabel('number of attempts','fontsize',15);
ylabel('number of jobs','fontsize',15);
fprintf('Type "return" to exit\n');
keyboard
