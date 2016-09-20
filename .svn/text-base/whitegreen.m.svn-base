function cm=whitegreen(gray_level)

if ~exist('gray_level','var')
  gray_level=0;
end

cm=flipud([ (1:(-gray_level/63):(1-gray_level))' repmat((0:(1-gray_level)/63:(1-gray_level))',1,2)]);
%cm=flipud([ (1-gray_level)*ones(64,1) repmat((0:(1-gray_level)/63:(1-gray_level))',1,2)]);
cm=cm(:,[2 1 3]);

if nargout==0
  colormap(cm);
end
