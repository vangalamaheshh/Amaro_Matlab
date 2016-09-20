function s=rand_uid(n)

s=[];
c=['0':'9' 'a':'z' 'A':'Z'];
for i=1:n
  s=[s c(ceil(rand*length(c)))];
end
s=char(s);
