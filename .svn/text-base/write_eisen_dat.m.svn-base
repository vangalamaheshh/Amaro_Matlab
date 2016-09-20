function write_eisen_dat(fname,gene_acc,gene_name,array_name,db_name, ...
			 dat,gweight,aweight,integer_val)

if exist('gene_acc','var') && iscell(gene_acc)
  gene_acc=strvcat(gene_acc);
end

if exist('gene_name','var') && iscell(gene_name)
  gene_name=strvcat(gene_name);
end

if exist('array_name','var') && iscell(array_name)
  array_name=strvcat(array_name);
end

if (nargin<9)
  integer_val=0;
end

if isempty(gweight)
  gw=0;
else 
  gw=1;
end

if isempty(aweight)
  aw=0;
else
  aw=1;
end

if integer_val
  dat=round(dat);
  dt='%d';
else
  dt='%f';
end

if strcmp(fname,'stdout')
  f=1;
else
  f=fopen(fname,'w');
end

fprintf(f,'%s\tNAME\t',db_name);
if gw 
  fprintf(f,'GWEIGHT\t');
end

for i=1:size(dat,2)
  if (i>1) 
    fprintf(f,'\t');
  end
  fprintf(f,'%s',deblank(array_name(i,:)));
end
fprintf(f,newline);

if aw
  fprintf(f,'EWEIGHT\t\t\t');
  for j=1:(length(aweight)-1)
    fprintf(f,'%f\t',aweight(j));
  end
  fprintf(f,['%f' newline],aweight(end));  
end

for i=1:size(dat,1)
  fprintf(f,'%s\t%s\t',deblank(gene_acc(i,:)),...
	  deblank(gene_name(i,:)));
  if gw
    fprintf(f,'%f\t',gweight(i));
  end
  for j=1:(size(dat,2)-1)
    if dat(i,j)==round(dat(i,j))
      fprintf(f,[ '%d' ' \t'],dat(i,j));
    else
      fprintf(f,[ dt ' \t'],dat(i,j));    
    end
  end
  if dat(i,end)==round(dat(i,end))
    fprintf(f,'%d',dat(i,end));
  else
    fprintf(f,dt,dat(i,end));
  end
  fprintf(f,newline);
end

if f~=1
  fclose(f);
end
