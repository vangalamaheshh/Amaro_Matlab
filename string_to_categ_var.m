function [categ_var]= string_to_categ_var(string_var)
categ_var=zeros(length(string_var),1);
[i m]=count(string_var);
for i=1:length(m)
    disp(sprintf('%d',i))
    categ_var(ismember(string_var,m{i}))=i;
end


end