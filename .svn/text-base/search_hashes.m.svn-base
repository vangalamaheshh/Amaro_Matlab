function [ord,s,res,h]=search_hashes(txt,hashes,aa4hash)

sw=0;
m=-1;
words=regexp(txt,'(\w*)','tokens');
res=cell(length(hashes),1);
for i=1:length(hashes)
  for k=1:length(words)
    w=lower(words{k}{1});
    tmp=get(hashes{i},w);
    if ~isempty(tmp) && tmp(1)~=-1
      res{i}{end+1,1}=w;
      res{i}{end,2}=tmp;
      m=max([tmp; m],[],1);
    elseif ~isempty(tmp) && tmp(1)==-1
      sw=1;
    end
  end
end

if (m>0)
  h=sparse(length(hashes),m);
  for i=1:length(res)
    for j=1:size(res{i},1)
      r=res{i}{j,2};
      tfidf=log2(m./length(r)); % actually only idf now
      h0=sparse(ones(length(r),1),r,tfidf*ones(length(r),1),1,m);
      h(i,:)=h(i,:)+h0;
    end
  end
  if exist('aa4hash','var')
    top_exact=~cellfun('isempty',regexp(lower(aa4hash), ...
                                        lower(regexprep(txt,'\(\S*\)','')))); 
    top_exact=sum(top_exact,2);
  else
    top_exact=[];
  end
  
  if isempty(top_exact)
    [o1,o2]=sort(-sum(h,1));
  else
    [o1,o2]=sort(-sum(h,1)-1e6*top_exact(1:size(h,2))');
  end    
  ord=full(o2(find(o1~=0)));
  s=-full(o1(find(o1~=0)));
  if (0) % old option
    if (sw==1) && exist('aa4hash','var') % we encountered stop words
      cur_hashes=make_hashes(aa4hash(ord,:));
      [ord2,s2,res2,h2]=search_hashes(txt,cur_hashes);
      ord=ord(ord2);
      %    s=s2;
      s=s(ord2)+0.1*s2;
      [ss,ssi]=sort(-s);
      s=s(ssi);
      ord=ord(ssi);
    end
  end
else
  h=[];
  ord=[];
  s=[];
end
