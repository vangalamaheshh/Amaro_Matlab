function [li,D,mv]=merge_levels(raw,rl,th)

md=0;
n=size(rl,1);
level=derunlength([rl(:,1:2) (1:n)']);
D=zeros(n,n);

for i=1:n
  li{i}=find(level==i);
end


% use upper triangle
for i=1:(n-1)
  for j=(i+1):n
    D(i,j)=level_dist(raw,li{i},li{j});
  end
end
D=D+tril(nan(n,n));

mv=[];
while md<th && ~isnan(md) && length(li)>1
  [md,mi]=min(D(:));
  mv(end+1)=md;
  disp(length(li));
  disp(md);
  if md<th && ~isnan(md)
    [c1,c2]=ind2sub(size(D),mi);
    if c1>=c2
      error('cant be');
    end
    
    % {c1 c2}-> c1
    disp(['merging [' num2str(median(raw(li{c1}))) ',' num2str(median(raw(li{c2}))) ']'  ...
          '->' num2str(median(raw([li{c1}; li{c2}])))]);
    li{c1}=[li{c1}; li{c2}];
    li(c2)=[];
    D(:,c2)=[];
    D(c2,:)=[];
    
    % update distances
    for i=1:(c1-1)
      D(i,c1)=level_dist(raw,li{i},li{c1});
    end
    for i=(c1+1):size(D,1)
      D(c1,i)=level_dist(raw,li{c1},li{i});
    end    
  end
end

function d=level_dist(raw,l1,l2)
%d=-log(ranksum(raw(l1),raw(l2))+eps);

m1=median(raw(l1));
m2=median(raw(l2));
s1=mad(raw(l1),1)*mad_factor;
s2=mad(raw(l2),1)*mad_factor;

%d=abs(m1-m2)./sqrt(s1^2/length(l1)+s2^2/length(l2)); % 20
d=abs(m1-m2)./(s1+s2);
