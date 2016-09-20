
function [a_s,ix]=shuffle_gct_arms(a)

ix = randperm(length(a));
for
a_s = a(ix,:);

end



