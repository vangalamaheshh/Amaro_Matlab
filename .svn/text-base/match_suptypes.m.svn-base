function [X,Y]=match_suptypes(X,xsi,Y,ysi)
X1=reorder_D_sup(X,'cols',xsi);
Y1=reorder_D_sup(Y,'cols',ysi);

[x_typeacc,x_typedesc,X2,x_range,x_non_empty]=decollapse_supdat(X1,xsi,-1);
X2=reorder_D_sup(X2,'cols',x_range); 
x_range=1:length(x_range);
[y_typeacc,y_typedesc,Y2,y_range,y_non_empty]=decollapse_supdat(Y1,ysi,-1);
Y2=reorder_D_sup(Y2,'cols',y_range);
y_range=1:length(y_range);

[Mt,m1,m2]=match_string_sets(Y2.supacc(y_range,:),X2.supacc(x_range,:));
extra_in_y=setdiff(y_range,m1);
extra_in_x=setdiff(x_range,m2);

if ~isempty(extra_in_y)
  X2=add_D_sup(X2,Y2.supacc(extra_in_y,:),Y2.supdesc(extra_in_y,:),...
               zeros(length(extra_in_y),size(X2.dat,2)),'cols');
end

if ~isempty(extra_in_x)
  Y2=add_D_sup(Y2,X2.supacc(extra_in_x,:),X2.supdesc(extra_in_x,:),...
               zeros(length(extra_in_x),size(Y2.dat,2)),'cols');
end

[Mt,m1,m2]=match_string_sets(Y2.supacc,X2.supacc);
assert((length(m1)==length(m2)) & (size(Y2.supacc,1)==size(X2.supacc,1)))

Y2=reorder_D_sup(Y2,'cols',m1);

X2=collapse_supdat(X2,x_typeacc,x_typedesc,1:size(X2.supacc,1),1);
Y2=collapse_supdat(Y2,y_typeacc,y_typedesc,1:size(Y2.supacc,1),1);

[X,sxi2]=add_D_sup(X,X2.supacc,X2.supdesc,X2.supdat,'cols');
X=reorder_D_sup(X,'col',setdiff(1:size(X.supdat,1),xsi));
[Y,syi2]=add_D_sup(Y,Y2.supacc,Y2.supdesc,Y2.supdat,'cols');
Y=reorder_D_sup(Y,'col',setdiff(1:size(Y.supdat,1),ysi));

% typeacc=strvcat(x_typeacc; 




