function [d_scans,no1,no2]=fit_scans(names,scans);

v1=zeros(size(names,1),1);
v2=zeros(length(scans),1);

for i=1:length(scans)
  f=findstrings(names,scans{i}.scan);
%  keyboard
  if ~isempty(f)
    if length(f)>1
      disp('confused')
    else
      d_scans{f(1)}=scans{i};
      v2(i)=1;
      v1(f(1))=1;
    end
  end
end

no1=find(v1==0);
no2=find(v2==0);
