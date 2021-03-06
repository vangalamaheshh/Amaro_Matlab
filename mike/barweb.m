function handles = barweb(barvalues, error_high, error_low, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)

master_linewidth = 1;

%
% Usage: handles = barweb(barvalues, errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type)
%
% Ex: handles = barweb(my_barvalues, my_errors, [], [], [], [], [], bone, [], bw_legend, 1, 'axis')
%
% barweb is the m-by-n matrix of barvalues to be plotted.
% barweb calls the MATLAB bar function and plots m groups of n bars using the width and bw_colormap parameters.
% If you want all the bars to be the same color, then set bw_colormap equal to the RBG matrix value ie. (bw_colormap = [1 0 0] for all red bars)
% barweb then calls the MATLAB errorbar function to draw barvalues with error bars of length error.
% groupnames is an m-length cellstr vector of groupnames (i.e. groupnames = {'group 1'; 'group 2'}).  For no groupnames, enter [] or {}
% The errors matrix is of the same form of the barvalues matrix, namely m group of n errors.
% Gridstatus is either 'x','xy', 'y', or 'none' for no grid.
% No legend will be shown if the legend paramter is not provided
% 'error_sides = 2' plots +/- std while 'error_sides = 1' plots just + std
% legend_type = 'axis' produces the legend along the x-axis while legend_type = 'plot' produces the standard legend.  See figure for more details
%
% The following default values are used if parameters are left out or skipped by using [].
% width = 1 (0 < width < 1; widths greater than 1 will produce overlapping bars)
% groupnames = '1', '2', ... number_of_groups
% bw_title, bw_xlabel, bw_ylabel = []
% bw_color_map = jet
% gridstatus = 'none'
% bw_legend = []
% error_sides = 2;
% legend_type = 'plot';
%
% A list of handles are returned so that the user can change the properties of the plot
% handles.ax: handle to current axis
% handles.bars: handle to bar plot
% handles.errors: a vector of handles to the error plots, with each handle corresponding to a column in the error matrix
% handles.legend: handle to legend
%
%
% See the MATLAB functions bar and errorbar for more information
%
% Author: Bolu Ajiboye
% Created: October 18, 2005 (ver 1.0)
% Updated: Dec 07, 2006 (ver 2.1)
% Updated: July 21, 2008 (ver 2.3)

% Get function arguments
minargs=3;
if nargin < minargs
	error('Must have at least the first two arguments:  barweb(barvalues, error_high error_low, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, barwebtype)');
elseif nargin == minargs
	width = 1;
	groupnames = 1:size(barvalues,1);
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+1
	groupnames = 1:size(barvalues,1);
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+2
	bw_title = [];
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+3
	bw_xlabel = [];
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+4
	bw_ylabel = [];
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+5
	bw_colormap = jet;
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+6
	gridstatus = 'none';
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+7
	bw_legend = [];
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+8
	error_sides = 2;
	legend_type = 'plot';
elseif nargin == minargs+9
	legend_type = 'plot';
end

change_axis = 0;
ymin = 0;
ymax = 0;

if any(size(error_high) ~= size(error_low)), error('error_high and error_low don''t match in size'); end
if size(barvalues,1) ~= size(error_high,1) || size(barvalues,2) ~= size(error_high,2)
	error('barvalues and errors matrix must be of same dimension');
else
	if size(barvalues,2) == 1
		barvalues = barvalues';
		error_high = error_high';
                error_low = error_low';
	end
	if size(barvalues,1) == 1
		barvalues = [barvalues; zeros(1,length(barvalues))];
		error_high = [error_high; zeros(1,size(barvalues,2))];
                error_low = [error_low; zeros(1,size(barvalues,2))];
		change_axis = 1;
	end
	numgroups = size(barvalues, 1); % number of groups
	numbars = size(barvalues, 2); % number of bars in a group
	if isempty(width)
		width = 1;
	end
	
	% Plot bars
	handles.bars = bar(barvalues, width,'edgecolor','k', 'linewidth', master_linewidth);
	hold on
	if ~isempty(bw_colormap)
		colormap(bw_colormap);
	else
		colormap(jet);
	end
	if ~isempty(bw_legend) && ~strcmp(legend_type, 'axis')
		handles.legend = legend(bw_legend, 'location', 'best', 'fontsize',12);
		legend boxoff;
	else
		handles.legend = [];
	end
	
	% Plot errors
	for i = 1:numbars
		x =get(get(handles.bars(i),'children'), 'xdata');
		x = mean(x([1 3],:));
		handles.errors(i) = errorbar(x, barvalues(:,i), error_high(:,i)-barvalues(:,i), error_low(:,i)-barvalues(:,i),...
                'k', 'linestyle', 'none', 'linewidth', master_linewidth);

                ymin = min(ymin,min([error_high(:,i);error_low(:,i)]));
		ymax = max(ymax,max([error_high(:,i);error_low(:,i)]));
	end
	
	if error_sides == 1
		set(gca,'children', flipud(get(gca,'children')));
	end
	
	ylim([ymin ymax*1.1]);
	xlim([0.5 numgroups-change_axis+0.5]);
	
	if strcmp(legend_type, 'axis')
		for i = 1:numbars
			xdata = get(handles.errors(i),'xdata');
			for j = 1:length(xdata)
				text(xdata(j),  -0.03*ymax*1.1, bw_legend(i), 'Rotation', 60, 'fontsize', 12, 'HorizontalAlignment', 'right');
			end
		end
		set(gca,'xaxislocation','top');
	end
	
	if ~isempty(bw_title)
		title(bw_title, 'fontsize',14);
	end
	if ~isempty(bw_xlabel)
		xlabel(bw_xlabel, 'fontsize',14);
	end
	if ~isempty(bw_ylabel)
		ylabel(bw_ylabel, 'fontsize',14);
	end
	
  set(gca, 'xticklabel', groupnames, 'box', 'off', 'ticklength', [0 0], 'fontsize', 12,...
     'xtick',1:numgroups, 'linewidth', master_linewidth,'xgrid','off','ygrid','off');

	if ~isempty(gridstatus) && any(gridstatus == 'x')
		set(gca,'xgrid','on');
	end
	if ~isempty(gridstatus) && any(gridstatus ==  'y')
		set(gca,'ygrid','on');
	end
	
	handles.ax = gca;
	
	hold off
end
