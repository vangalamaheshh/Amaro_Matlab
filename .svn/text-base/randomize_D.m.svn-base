function Dr=randomize_D(D,opt)

Dr=D;
Dr.random=1;
switch opt
 case 'cols in each row'
  disp(1);
  for i=1:size(D.dat,1)
    r=randperm(size(D.dat,2));
    Dr.dat(i,:)=D.dat(i,r);
  end
 case 'cols'
  disp(2);
  r=randperm(size(D.dat,2));
  Dr.dat=D.dat(:,r);
end

