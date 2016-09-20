function s = fullmin(m)
if length(m)==1, s = m;
else s = fullmin(min(m)); end
end
