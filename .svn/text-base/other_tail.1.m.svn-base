function res=other_tail(pv,ipv)

if(0)
fwd=1;
back=length(pv);
%ipv=1-pv;

res=zeros(length(pv),1);
while(fwd<back)
  disp(['f=' num2str(fwd) '; b=' num2str(back)]);
  if pv(fwd)<(ipv(back))
    res(back)=pv(fwd);
    fwd=fwd+1;
  elseif pv(fwd)>(ipv(back))
    res(fwd)=ipv(back);
    back=back-1;
  else
    res(fwd)=ipv(back);
    res(back)=pv(fwd);
    fwd=fwd+1;
    back=back-1;
  end
end
disp(['f=' num2str(fwd) '; b=' num2str(back)]);
else
  res=zeros(length(pv),1);
  ipv=1-pv;
  for i=1:length(pv)
    if pv(i)<=0.5
      res(i)=ipv(min(find(ipv<=pv(i)+eps)));
    else
      res(i)=pv(max(find(pv<=ipv(i)+eps)));
    end
  end
end




