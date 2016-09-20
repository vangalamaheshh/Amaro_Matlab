function xy=getxy(m)

figure(1); clf
axis([0 1 0 1]);
hold on;
n=1;
while m>0
  for i=1:m 
    [xi,yi] = ginput(1);
    disp([n xi yi]);
    plot(xi,yi,'ro'); hold on;
    n = n+1;
    xy(:,n) = [xi;yi];  
  end
  m=input('more points?');
end

