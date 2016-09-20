function c=num2affy_call(v)
b=v;
b(v==0)='A';
b(v==1)='M';
b(v==2)='P';
b=char(b);
c=mat2cell(b,ones(length(v),1),1);

return

c=cell(size(s,1),1);
for i=1:length(c)
  c{i}=s(i,1);
end


s=num2str(v);
s(s=='0')='A';
s(s=='1')='M';
s(s=='2')='P';

c=cell(size(s,1),1);
for i=1:length(c)
  c{i}=s(i,1);
end

