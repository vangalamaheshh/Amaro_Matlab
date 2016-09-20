function X=remove_unknown_segments(X,len,only_ret_unknown_ret)

if exist('only_ret_unknown_ret','var')
  only_ret_unknown_ret=0;
end

dat=X.dat;
dat(isnan(dat))=0;
dif=diff([-ones(1,size(dat,2)); dat],1);
dat2=dat;
for i=1:size(dat,2)
  st=find(dif(:,i)~=0);
  
  if length(st)>2
    en=[ st(2:end)-1; size(dat,1) ];
    
    st=st(2:(end-1));
    en=en(2:(end-1));
    
    % Take care of beginning and ends of chromosomes.
    % maybe, allow len to be in Mb units
    
    s131=find((en-st+1<=len) & ( (dat(st-1,i)==1) & (dat(st,i)==3) & ...
                                 (dat(en+1,i)==1)));
    s232=find((en-st+1<=len) & ( (dat(st-1,i)==2) & (dat(st,i)==3) & ...
                                 (dat(en+1,i)==2)));
    
    for j=1:length(s131)
      dat2(st(s131(j)):en(s131(j)),i)=1;
    end
    
    if ~only_ret_unknown_ret
      for j=1:length(s232)
        dat2(st(s232(j)):en(s232(j)),i)=2;
      end
    end
    
  end
end

X.smooth=dat2;

