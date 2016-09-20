function [rdesc,racc,rweight,cdesc,cweight,dbname,dat]=read_eisen_dat(fname)

l=my_textread(fname,'%s','whitespace','\n','bufsize',1000000);
pos=findstr(9,l{1});
gw=0; aw=0;
rweight=[];
cweight=[];

dbname=l{1}(1:(pos(1)-1));
if (strcmp(l{1}(pos(2)+1:pos(3)-1),'GWEIGHT'))
   gw=1;
end
cdesc=[];
nc=size(pos,2)-1;
for j=3:(nc+1)
    cur_cdesc=l{1}((pos(j-1)+1):(pos(j)-1));
    if isempty(cur_cdesc)
        cdesc=strvcat(cdesc,'EMPTY');
    else
        cdesc=strvcat(cdesc,cur_cdesc);
    end
end
cdesc=strvcat(cdesc,l{1}((pos(j)+1):end));

pos=findstr(9,l{2});
if strcmp(l{2}(1:pos(1)-1),'EWEIGHT')
   aw=1;
end

rdesc=[];
racc=[];
nr=size(l,1)-1;
dat=zeros(nr,nc);
for k=2:size(l,1)
   if mod(k,100)==0
      disp(k)
   end
   pos=findstr(9,l{k});
   racc=strvcat(racc,l{k}(1:pos(1)-1));
   cur_rdesc=l{k}((pos(1)+1):(pos(2)-1));
   if isempty(cur_rdesc)
        rdesc=strvcat(rdesc,'EMPTY');
   else
        rdesc=strvcat(rdesc,cur_rdesc);
   end
   for j=1:nc
      if (j<nc)
	s=l{k}((pos(j+1)+1):(pos(j+2)-1));
      else
	s=l{k}((pos(j+1)+1):end);
      end
      if isempty(s) || strcmp(lower(deblank(s)),'na') || strcmp(lower(deblank(s)),'nan') || strcmp(lower(deblank(s)),'#n/a') 
         dat(k-1,j)=NaN;
      else
	 try
           dat(k-1,j)=sscanf(s,'%f');
	 catch
	   disp([ lasterr '[ ' s  ']@(' num2str(k-1) ',' num2str(j) ')'] );
	 end
      end
   end
end   
if gw
  if aw 
    cweight=dat(1,2:end);
  else
    cweight=dat(1,:);
  end    
  cdesc=cdesc(2:end,:);
  dat=dat(2:end,:);
end
if aw
  rweight=dat(:,1);
  rdesc=rdesc(2:end,:);
  racc=racc(2:end,:);
  dat=dat(:,2:end);
end

if nargout==1
  D.dat=dat;
  D.gacc=racc;
  D.gdesc=rdesc;
  D.sdesc=cdesc;
  D.sweight=cweight;
  D.gweight=rweight;
  D.dbname=dbname;
  rdesc=D;
end


