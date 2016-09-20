function s=disp_confmat(D,confmat,suprange)

n=sum(confmat,1);
ncorr=diag(confmat);
nerr=n-ncorr';
s=[ sprintf(['Type\t#Samples\t#Errors\tBreakdown of errors' newline])];
for i=1:size(confmat,1)
  s=[ s sprintf('%s\t%d\t%d',deblank(D.supacc(suprange(i),:)),n(i),nerr(i))];
  st=[];
  if nerr(i)>0
    x=confmat(:,i);
    x(i)=0;
    xpos=find(x>0);
    [sx,si]=sort(x(xpos),1,'descend');
    for j=1:length(si) 
      st=[ st sprintf('%d %s,',sx(j),deblank(D.supacc(suprange(xpos(si(j))),:)))];
    end
    st=st(1:(end-1));
  end
  s=[ s sprintf(['\t%s' newline],st)];
end

if nargout==0
  disp(s);
end
