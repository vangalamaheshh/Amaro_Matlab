function roc_curve(func_handle,data,cls,pranges)

xi=zeros(length(pranges),1);

while( xi(end)<length(pranges{end}) )
  r=1;
  for i=1:length(xi)
    xi(1)=xi(1)+r;
    if (xi(i) >= length(pranges{i}))
      xi(i)=0;
      r=1;
    else
      r=0;
    end
    disp(xi);
  end
  for i=1:length(xi)
    x(i)=pranges{i}(xi(i)+1);
  end
  [fp,fn,feval(func_handle(data,cls,x));
end

