function comp_dist(x,nstep,is_dist)
% comp_dist(x,nstep,is_dist)

if nargin == 2
  is_dist = 0;
end

n=length(x);
xmin=[];
xmax=[];
if length(nstep) == 1
  for i=1:n
    xmin=min([xmin; x{i}]);
    xmax=max([xmax; x{i}]);
  end
  step=(xmax-xmin)/nstep;
  xpos=xmin:step:xmax;
else
  xpos=nstep;
end

if length(xpos)<1
    bar(1,1);
    text(1,1.3,'Problem with figure');
    axis([0 2 0 2]);
    return
end

hs=zeros(n,length(xpos));
for i=1:n
  [h,a]=hist(x{i},xpos);
  hs(i,:)=h;
  % in the future might use histc
%   h=histc(x{i},xpos);
%   hs(i,:)=h;
end
if is_dist>0
  sums=sum(hs,2);
  hs = hs ./ repmat(sums,1,length(xpos));
  if is_dist==2
    hs=cumsum(hs,2);
  end
end

bar(xpos,hs',1,'grouped');
colormap(cool)
if is_dist
  ylabel('probability')
  legend([ num2str((1:n)') repmat(' n=',n,1) num2str(sums) ])
else
  ylabel('count');
  legend(num2str((1:n)') )
end

