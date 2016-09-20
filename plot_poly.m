function ph=plot_poly(pol,n_points,varargin)

if ~exist('n_points','var') || isempty(n_points)
  n_points=1000;
end

ax=axis;
x=ax(1):(ax(2)-ax(1))/(n_points-1):ax(2);

pol=fliplr(pol);
X=zeros(length(pol),size(x,2));
for i=1:length(pol)
  X(i,:)=x.^(i-1);
end
y=pol*X;

ph=plot(x,y,varargin{:});
axis(ax);

