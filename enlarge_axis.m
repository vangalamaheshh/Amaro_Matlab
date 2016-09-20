function enlarge_axis(p)
if length(p)==1
  p=[p p p; p p p];
end
if size(p,1)==1
  p=[p; p];
end

ax=axis;
Xw=ax(2)-ax(1);
Yw=ax(4)-ax(3);
if length(ax)>4
  Zw=ax(6)-ax(5);
  nax=[ ax(1)-Xw*p(1,1) ax(2)+Xw*p(2,1) ax(3)-Yw*p(1,2) ax(4)+Yw*p(2,2) ax(5)-Zw*p(1,3) ax(6)+Zw*p(2,3)];
else
  nax=[ ax(1)-Xw*p(1,1) ax(2)+Xw*p(2,1) ax(3)-Yw*p(1,2) ax(4)+Yw*p(2,2)];
end
axis(nax);

