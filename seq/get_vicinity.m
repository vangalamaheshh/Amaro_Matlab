function v = get_vicinity(chr,pos,radius)

if ~exist('radius','var'), radius=5; end

v = cell(length(chr),1);
for i=1:length(chr), if ~mod(i,100), fprintf('%d/%d ',i,length(chr)); end
  v{i} = genome_region(chr(i),pos(i)-radius,pos(i)+radius);
end,fprintf('\n');

