function bar_with_colors(varargin)
% bar_with_colors(varargin)
%
% all parameters but the last one are passed unchanged to bar()
% last parameter should be matrix C, the same x*y size as the matrix being plotted, with size(C,3)==3,
% with page1,2,3=RGB for each matrix element
%
% Mike Lawrence 2011-01-27


fprintf('WARNING: DOES NOT WORK PROPERLY: colors are out of synch\n'); %keyboard;

h = bar(varargin{1:end-1});


ch = get(h,'children');
fvd = get(ch,'Faces');
fvcd = get(ch,'FaceVertexCData');

C = varargin{end};

colormap(C);

z = get(ch,'FaceVertexCData');
if length(z)==size(C,1)
  set(ch,'FaceVertexCData',(1:size(C,1))');
elseif length(z)==size(C,1)*5 + 1
  q = round(0:0.2:size(C,1))';
  set(ch,'FaceVertexCData',q);
else
  fprintf('Help!\n');
end

% (hacky implementation, might not work in all cases)
