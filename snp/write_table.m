function write_table(fname,tbl,sep,empty_element)

if ~exist('sep','var') || isempty(sep)
  sep={char(9),sprintf(newline)};
end
if ~exist('empty_element','var')
  empty_element='---';
end

c=1;
r=1;

if ischar(fname)
  f=fopen(fname,'w');
elseif f>0
  f=fname;
end


for r=1:size(tbl,1)
  for c=1:size(tbl,2)
    x=tbl{r,c};
    if isempty(x)
      fprintf(f,'%s',empty_element);
    elseif ischar(x)
      fprintf(f,'%s',x);
    elseif isnumeric(x)
      if prod(size(x))==1
        if x==round(x)
          fprintf(f,'%d',x);
        else
          fprintf(f,'%f',x);
        end
      else
        x=mat2cell(x,ones(1,size(x,1)),ones(1,size(x,2)));
        write_table(f,x,sep(2:end,:));
      end
    elseif iscell(x)
      write_table(f,x,sep(2:end,:));
    else
      error('unexpected type');
    end
    if c<size(tbl,2)
      fprintf(f,'%s',sep{1,1});
    end
  end
  if r<size(tbl,1)
    fprintf(f,'%s',sep{1,2});
  end
end

if ischar(fname)
  fprintf(f,'%s',sep{1,2});
  fclose(f);
end
