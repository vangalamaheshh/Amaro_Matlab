function ranks = findranks(input,dim,direction);

if ~exist('dim','var'), dim=1; end
if ~exist('direction','var'), direction='ascend'; end

[tmp ord] = sort(input,dim,direction);
[tmp ranks] = sort(ord,dim,'ascend');
