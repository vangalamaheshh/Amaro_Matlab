function p=draw_block(x,y,z,c)

d=0.7;
p(1)=patch([x(1) x(2) x(2) x(1) x(1)],[ y(1) y(1) y(2) y(2) y(1)],repmat(z(1),1,5),c*d);
p(2)=patch([x(1) x(2) x(2) x(1) x(1)],[ y(1) y(1) y(2) y(2) y(1)],repmat(z(2),1,5),c*d);
p(3)=patch(repmat(x(1),1,5),[ y(1) y(1) y(2) y(2) y(1)],[z(1) z(2) z(2) z(1) z(1)],c);
p(4)=patch(repmat(x(2),1,5),[ y(1) y(1) y(2) y(2) y(1)],[z(1) z(2) z(2) z(1) z(1)],c);
p(5)=patch([x(1) x(2) x(2) x(1) x(1)],repmat(y(1),1,5),[z(1) z(1) z(2) z(2) z(1)],d*c+(1-d)*[1 1 1]);
p(6)=patch([x(1) x(2) x(2) x(1) x(1)],repmat(y(2),1,5),[z(1) z(1) z(2) z(2) z(1)],d*c+(1-d)*[1 1 1]);

