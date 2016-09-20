function legos(varargin)
% legos(input1,input2,input3)
% legos({input1,input2,input3},titles)
% legos({input1,input2,input3},titles,P)

if nargin==1 && iscell(varargin{1})
  z = varargin{1};
  titles = [];
  P=[];
elseif nargin==2 && iscell(varargin{1}) && iscell(varargin{2})
  z = varargin{1};
  titles = varargin{2};
  P=[];
elseif nargin==3 && iscell(varargin{1}) && iscell(varargin{2}) && isstruct(varargin{3})
  z = varargin{1};
  titles = varargin{2};
  P = varargin{3};
else
  z = varargin;
  titles = [];
  P=[];
end

P = impose_default_value(P,'nrows',[]);

n = length(z);
if ~isempty(P.nrows)
  nrows=P.nrows;
  ncols = ceil(n/nrows);
else
  if n<=3
    nrows=1;
    ncols=n;
  else
    s = ceil(sqrt(n));
    nrows=s;
    ncols=s;
  end
end

fi=1;
for y=1:nrows, for x=1:ncols
    if fi<=n
      subplot('position',[(x-1)*(1/ncols) 1-(y*(1/nrows)) (1/ncols) (1/nrows)]);
      lego(z{fi},P);
      if ~isempty(titles)
        zl=zlim;
        text(1,1,zl(2)*1.6,titles{fi},'fontsize',8);
      end
    end
    fi=fi+1;
end,end,set(gcf,'color',[1 1 1]);

