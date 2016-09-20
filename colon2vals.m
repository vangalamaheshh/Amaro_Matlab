function vals=colon2vals(st)

pos=find(st==':');
lpos=length(pos);
if lpos==0
  vals=str2num(st);
elseif lpos==1
  vals=str2num(st(1:(pos(1)-1))):str2num(st((pos(1)+1):end));
elseif lpos==2
  vals=str2num(st(1:(pos(1)-1))):str2num(st((pos(1)+1):(pos(2)-1))):str2num(st((pos(2)+1):end));
else
  error('too many :');
end
  
