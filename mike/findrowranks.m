function ranks = findrowranks(input,dir)
if ~exist('dir','var'), dir = 'ascend'; end
ranks = findranks(input,2,dir);
