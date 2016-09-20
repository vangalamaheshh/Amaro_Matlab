
function [a_s,ix]=shuffle_gct_arms(a)
a_s=zeros(size(a));
for i =1:size(a,2)
    ix = randperm(size(a,1));
    a_s(:,i) = a(ix,i);
end

end



