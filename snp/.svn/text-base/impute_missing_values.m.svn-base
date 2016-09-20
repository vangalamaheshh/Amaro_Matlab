function C1=impute_missing_values(C,params)

if ischar(params)
  tmp.method=params;
  params=tmp;
end

C1=C;
switch params.method
 case 'zero'
  C1.dat(isnan(C1.dat))=0;
 case 'before'
  C1=C;
  rl=runlength(isnan(C.dat));
  if ~iscell(rl)
    rl={rl};
  end
  for j=1:length(rl)
    seg=find(rl{j}(:,3)==1);
    if ~isempty(seg)
      for k=1:length(seg)
        if rl{j}(seg(k),1)==1
          C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=0;
        else
          C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=C.dat(rl{j}(seg(k),1)-1,j);
        end
      end
    end
  end
  
 case 'ifsame'
  C1=C;
  rl=runlength(isnan(C.dat),C.chrn);
  if ~iscell(rl)
    rl={rl};
  end
  for j=1:length(rl)
    disp(j)
    seg=find(rl{j}(:,3)==1);
    if ~isempty(seg)
      for k=1:length(seg)
        if rl{j}(seg(k),1)~=1 && rl{j}(seg(k),2)~=size(C.dat,1) && ...
              C.chrn(rl{j}(seg(k),1)-1)==C.chrn(rl{j}(seg(k),1)) && ...
              C.chrn(rl{j}(seg(k),2)+1)==C.chrn(rl{j}(seg(k),2)) && ...
              C.dat(rl{j}(seg(k),1)-1,j)==C.dat(rl{j}(seg(k),2)+1,j)
          C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=C.dat(rl{j}(seg(k),1)-1,j);
        end
      end
    end
  end
  
 case 'max'
  C1=C;
  for j=1:size(C1.dat,2)
    if nnz(isnan(C.dat(:,j)))>0
      rl{j}=runlength(isnan(C.dat(:,j)),C.chrn);
      seg=find(rl{j}(:,3)==1);
      if ~isempty(seg)
        for k=1:length(seg)
          if rl{j}(seg(k),1)==1 || C.chrn(rl{j}(seg(k),1)-1)~=C.chrn(rl{j}(seg(k),1))
            if rl{j}(seg(k),2)==size(C.dat,1) || C.chrn(rl{j}(seg(k),2)+1)~=C.chrn(rl{j}(seg(k),2))
              C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=0;
            else
              C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=C.dat(rl{j}(seg(k),2)+1,j);
            end
          elseif rl{j}(seg(k),2)==size(C.dat,1) || C.chrn(rl{j}(seg(k),2)+1)~=C.chrn(rl{j}(seg(k),2))
            if rl{j}(seg(k),1)==1 || C.chrn(rl{j}(seg(k),1)-1)~=C.chrn(rl{j}(seg(k),1))
              C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=0;
            else
              C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=C.dat(rl{j}(seg(k),1)-1,j);
            end
          else
            C1.dat(rl{j}(seg(k),1):rl{j}(seg(k),2),j)=max(C.dat(rl{j}(seg(k),1)-1,j),C.dat(rl{j}(seg(k),2)+1,j));
          end
        end
      end
    end
  end
 

 case 'interp'
  for i=1:max(C.chrn)
    disp(i);
    in_chr=find(C.chrn==i);
    rl=runlength(isnan(C.dat(in_chr,:)));
    for j=1:length(rl)
      seg=find(rl{j}(:,3)==1);
      for k=1:length(seg)
        if rl{j}(seg(k),1)==1 % beginning of chr
          if rl{j}(seg(k),2)<length(in_chr)
            C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=C1.dat(in_chr(1)-1+rl{j}(seg(k),2)+1,j);
          else
            C1.dat(in_chr,j)==0;
          end
        elseif rl{j}(seg(k),2)==length(in_chr) % end of chr
          C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=C1.dat(in_chr(1)-1+rl{j}(seg(k),1)-1,j);        
        else
          a=C1.dat(in_chr(1)-1+rl{j}(seg(k),1)-1,j);
          b=C1.dat(in_chr(1)-1+rl{j}(seg(k),2)+1,j);
          x=C1.pos(in_chr(1)-1+((rl{j}(seg(k),1)-1):rl{j}(seg(k),2)+1));
          y=(b*(x-x(1))+a*(x(end)-x))./(x(end)-x(1));
          C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=y(2:(end-1));
        end
      end
    end
  end
 
 case 'maxabs'
  for i=1:max(C.chrn)
    disp(i);
    in_chr=find(C.chrn==i);
    rl=runlength(isnan(C.dat(in_chr,:)));
    for j=1:length(rl)
      seg=find(rl{j}(:,3)==1);
      for k=1:length(seg)
        if rl{j}(seg(k),1)==1 % beginning of chr
          if rl{j}(seg(k),2)<length(in_chr)
            C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=C1.dat(in_chr(1)-1+rl{j}(seg(k),2)+1,j);
          else
            C1.dat(in_chr,j)==0;
          end
        elseif rl{j}(seg(k),2)==length(in_chr) % end of chr
          C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=C1.dat(in_chr(1)-1+rl{j}(seg(k),1)-1,j);        
        else
          a=C1.dat(in_chr(1)-1+rl{j}(seg(k),1)-1,j);
          b=C1.dat(in_chr(1)-1+rl{j}(seg(k),2)+1,j);
          if abs(a)>abs(b)
            C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=a;
          else
            C1.dat(in_chr(1)-1+(rl{j}(seg(k),1):rl{j}(seg(k),2)),j)=b;
          end
        end
      end
    end
  end
end
