function [hashes,stopwords]=make_hashes(lists)

hashes=cell(size(lists,2),1);
stopwords={};
for i=1:size(lists,2)
  hashes{i}=java.util.Hashtable;
  for j=1:size(lists,1)
    if mod(j,1000)==0
      disp(j);
    end
    words=regexp(lists{j,i},'(\w*)','tokens');
    for k=1:length(words)
      w=lower(words{k}{1});
      cur=get(hashes{i},w);
      if length(cur)>1000 % && isempty(regexp(w,'(\d*)'))
          put(hashes{i},w,-1);
          stopwords{end+1}=w;
          stopwords
      elseif isempty(cur)
        put(hashes{i},w,j);
      elseif cur(1)~=-1
        put(hashes{i},w,[cur; j]);
      end
    end
  end
end 

