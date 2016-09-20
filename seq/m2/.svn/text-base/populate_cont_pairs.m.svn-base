function [ x ] = populate_cont_pairs(cov_track, cat_track, cont_new)
% cont_new = context-newbase
%
% returns x:  each row = [position newbase] of a throwable mutation

x = cell(size(cont_new, 1), 1);

for i = 1:size(cont_new,1)
    idx = find(cat_track == cont_new(i, 1));
    y = cell(length(idx), 1);
    for j = 1:length(idx)
        y{j} = repmat([idx(j), cont_new(i, 2)], cov_track(idx(j)), 1);
    end
    x{i} = cat(1, y{:});
end

x = cat(1, x{:});


