function test_colormaps

nms={'GRAY','HOT', 'COOL', 'BONE', 'COPPER', 'PINK', 'FLAG', 'PRISM', 'JET'};

for i=1:length(nms)
  colormap(lower(nms{i}));
  disp(lower(nms{i}));
  pause;
end
