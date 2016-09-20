function r=count_in_base(dig,base)

N=base^dig;
n1=1;
n2=N/base;
r=[];
for i=1:dig
   tmp=repmat(1:base,n1,n2);
   r=[ tmp(:) r];
   n1=n1*base;
   n2=n2/base;
end


