function distribution_ordered_plot(varargin)
figure()
hold on
for i=1:length(varargin)
    step=.2/(length(varargin{i})-1);
    
    scatter((i-.1):step:((i+.1)),sort(varargin{i}),40,[0 0 0],'filled')
    line([(i-.1):.05:(i+.1)],repmat(mean(varargin{i}(~isnan(varargin{i}))),5),'Color',[1 0 0]);

end
end