function s = fullmean(m)
if length(m)==1, s = m;
else s = fullmean(mean(m)); end
end
