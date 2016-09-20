function Y = collapse_categories_to_65_tr(X,categs_txt)
% collapses rows (by summing) to the standard context65 list,
% using the specified categs_txt as a guide

if ndims(X)>2, error('doesn''t work with 3D matrices'); end

k = map_categories_to_65(categs_txt);

if any(isnan(k)) || any(k<1) || any(k>65) || any(k~=round(k))
  error('problem with map_categories_to_65');
end
if length(k)~=size(X,2), error('X has %d rows, whereas %s has %d categories',size(X,2),categs_txt,length(k)); end

Y = zeros(65,size(X,1));
for i=1:length(k)
  Y(k(i),:) = Y(k(i),:) + X(:,i)';
end
