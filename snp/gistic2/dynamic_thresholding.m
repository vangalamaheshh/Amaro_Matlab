verbose('Setting thresholds dynamically...',30);
QA = Qs.amp;
QD = Qs.del;
QD(:,4) = -1*QD(:,4);
amplitudes = [QA(:,4); QD(:,4)];

hx = -2:0.01:2;
hc = histc(amplitudes,hx);
qq = smooth(hc/sum(hc));
[a,b,m,mi] = extrema(qq);
right_lim = min(mi(find(mi > median(1:length(qq)))));
left_lim = max(mi(find(mi < median(1:length(qq)))));

t_amp = hx(right_lim);
params.t_amp = t_amp;
t_del = hx(left_lim);
params.t_del = t_del;

verbose(['t_amp = ' num2str(params.t_amp)],30);
verbose(['t_del = ' num2str(params.t_del)],30);

figure()
subplot(2,5,[1:3 6:8])
bar(hx,hc/sum(hc))
xlim([min(hx) max(hx)])
xlabel('Difference in segmentation levels')
ylabel('Frequency')
line([params.t_amp params.t_amp],[0 max(ylim)],'Color','green','LineStyle','--');
line([params.t_del params.t_del],[0 max(ylim)],'Color','green','LineStyle','--');

subplot(2,5,4:5)
bar(hx(right_lim-15:right_lim+15),hc(right_lim-15:right_lim+15)/ ...
    sum(hc))
xlabel('Difference in segmentation levels')
ylabel('Frequency')
title(['Amp boundary = ' num2str(params.t_amp)])
xlim([hx(right_lim-15) hx(right_lim+15)])
line([params.t_amp params.t_amp],[0 max(ylim)],'Color','green','Linestyle','--');

subplot(2,5,9:10)
bar(hx(left_lim-15:left_lim+15),hc(left_lim-15:left_lim+15)/ ...
    sum(hc))
xlabel('Difference in segmentation levels')
ylabel('Frequency')
title(['Del boundary = ' num2str(params.t_del)])
xlim([hx(left_lim-15) hx(left_lim+15)])
line([params.t_del params.t_del],[0 max(ylim)],'Color','green','Linestyle','--');
    
saveas(gcf,[base_dir 'dynamic_threshold_histogram.fig']);
saveas(gcf,[base_dir 'dynamic_threshold_histogram.pdf'],'pdf');
t_del = -1*t_del;
params.t_del = -1*params.t_del;
