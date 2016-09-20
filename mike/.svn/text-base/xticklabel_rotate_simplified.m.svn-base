function hText = xticklabel_rotate_simplified(XTick,rot,varargin)
%xticklabel_rotate_simplified(XTick,rot,varargin)
%
% based on the function below, originally by Brian FG Katz
%
% modified by Mike Lawrence 2009-07-29
%   --> Removed the latter half of the procedure,
%         which did automatic resizing of the plot area.
%       This was troublesome when using "matlab -nodisplay",
%         because under those conditions, get(gca,'position')
%         returns negative values for height and/or width.
%         (is this a MATLAB bug??)
%
%hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
%
% Syntax: xticklabel_rotate
%
% Input:    
% {opt}     XTick       - vector array of XTick positions & values (numeric) 
%                           uses current XTick values or XTickLabel cell array by
%                           default (if empty) 
% {opt}     rot         - angle of rotation in degrees, 90by default
% {opt}     XTickLabel  - cell array of label strings
% {opt}     [var]       - "Property-value" pairs passed to text generator
%                           ex: 'interpreter','none'
%                               'Color','m','Fontweight','bold'
%
% Output:   hText       - handle vector to text labels
%
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45and change
% font size
%    xticklabel_rotate([],45,[],'Fontsize',14)
%
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Note : you can not re-run xticklabel_rotate on the same graph. 
%
% 


% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%                       Arbitrary angle rotation
%                       Output of text handles
%                       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size 
%                           and keep text on plot
%                           (handles small window resizing after, but not well due to proportional placement with 
%                           fixed font size. To fix this would require a serious resize function)
%                       Uses current XTick by default
%                       Uses current XTickLabel is different from XTick values (meaning has been already defined)

% Brian FG Katz
% bfgkatz@hotmail.com
% 23-05-03
% Modified 03-11-06 after user comment
%Allow for exisiting XTickLabel cell array

% Other m-files required: cell2mat
% Subfunctions: none
% MAT-files required: none
%
% See also: xticklabel_rotate90, TEXT,  SET

% Based on xticklabel_rotate90
%   Author: Denis Gilbert, Ph.D., physical oceanography
%   Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%   email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%   February 1998; Last revision: 24-Mar-2003

% check to see if xticklabel_rotate has already been here (no other reason for this to happen)
if isempty(get(gca,'XTickLabel')),
    error('xticklabel_rotate : can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased')  ;
end

% if no XTickLabel AND no XTick are defined use the current XTickLabel
%if nargin < 3 & (~exist('XTick') | isempty(XTick)),
if (nargin < 3 || isempty(varargin{1})) & (~exist('XTick') | isempty(XTick)),
  xTickLabels = get(gca,'XTickLabel')  ; % use current XTickLabel
                                         if ~iscell(xTickLabels)
                                           % remove trailing spaces if exist (typical with auto generated XTickLabel)
                                           temp1 = num2cell(xTickLabels,2)         ;
                                           for loop = 1:length(temp1),
                                             temp1{loop} = deblank(temp1{loop})  ;
                                             end
                                             xTickLabels = temp1                     ;
                                             end
                                             
end

% if no XTick is defined use the current XTick
if (~exist('XTick') | isempty(XTick)),
    XTick = get(gca,'XTick')        ; % use current XTick 
end

%Make XTick a column vector
XTick = XTick(:);

if ~exist('xTickLabels'),
  % Define the xtickLabels 
  % If XtickLabel is passed as a cell array then use the text
  if (length(varargin)>0) & (iscell(varargin{1})),
        xTickLabels = varargin{1};
        varargin = varargin(2:length(varargin));
        else
        xTickLabels = num2str(XTick);
        end
end    

if length(XTick) ~= length(xTickLabels),
    error('xticklabel_rotate : must have same number of elements in "XTick" and "XTickLabel"')  ;
end

%Set the Xtick locations and set XTicklabel to an empty string
set(gca,'XTick',XTick,'XTickLabel','')

if nargin < 2,
    rot = 90 ;
end

% Determine the location of the labels based on the position
% of the xlabel
hxLabel = get(gca,'XLabel');  % Handle to xlabel
xLabelString = get(hxLabel,'String');

% if ~isempty(xLabelString)
%    warning('You may need to manually reset the XLABEL vertical position')
% end

set(hxLabel,'Units','data');
xLabelPosition = get(hxLabel,'Position');
y = xLabelPosition(2);

%CODE below was modified following suggestions from Urs Schwarz
y=repmat(y,size(XTick,1),1);
% retrieve current axis' fontsize
fs = get(gca,'fontsize');

% Place the new xTickLabels by creating TEXT objects
hText = text(XTick, y, xTickLabels,'fontsize',fs);

% Rotate the text objects by ROT degrees
set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
