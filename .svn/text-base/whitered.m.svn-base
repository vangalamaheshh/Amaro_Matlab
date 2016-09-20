function cm=whitered(gray_level)

if ~exist('gray_level','var')
  gray_level=eps;
end

cm=flipud([ (1:(-gray_level/63):(1-gray_level))' repmat((0:(1-gray_level)/63:(1-gray_level))',1,2)]);
%cm=flipud([ (1-gray_level)*ones(64,1) repmat((0:(1-gray_level)/63:(1-gray_level))',1,2)]);

if nargout==0
  colormap(cm);
end
