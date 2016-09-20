function distribution_ordered_plot_deTiN(varargin)
figure()
hold on
for i=1:length(varargin)
    step=.2/(length(varargin{i})-1);
    sum(varargin{i}<.02);
    x=[i-.1:step:i+.1];
    [Y I]=sort(varargin{i});
        scatter(x(Y>=.02),Y(Y>=.02),40,[1 0 0],'filled')

    scatter(x(Y<.02),Y(Y<.02),40,[0 0 1],'filled')
    line([(i-.1):.05:(i+.1)],repmat(mean(varargin{i}(~isnan(varargin{i}))),5),'Color',[0 0 0]);

end
end