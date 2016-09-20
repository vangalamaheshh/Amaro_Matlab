function write_cdt_gtr_atr(fname,ofname,lnk,idx,clnk,cidx)
N=length(idx);
if isempty(lnk)
  lnk=[ ones(N-1,1) (1:(N-1))' (2:N)' (2:N)' ones(N-1,1)];
end

nd=zeros(N,1);
f=fopen([fname '.gtr'],'w');
for i=1:size(lnk,1)
   fprintf(f,'NODE%dX\t',i);
   if nd(lnk(i,1))
      fprintf(f,'NODE%dX\t',nd(lnk(i,1)));
   else
      fprintf(f,'GENE%dX\t',idx(lnk(i,1)));
   end
   if nd(lnk(i,3))
      fprintf(f,'NODE%dX\t',nd(lnk(i,3)));
   else
      fprintf(f,'GENE%dX\t',idx(lnk(i,3)));
   end
   nd(lnk(i,1))=i;
   fprintf(f,'%f%s',(lnk(i,5)/(max(lnk(:,5)+1))*1.8-0.8),char([13 10]));
end
fclose(f);

cN=length(cidx);
if isempty(clnk)
  clnk=[ ones(cN-1,1) (1:(cN-1))' (2:cN)' (2:cN)' ones(cN-1,1)];
end
cnd=zeros(cN,1);
f=fopen([fname '.atr'],'w');
for i=1:size(clnk,1)
   fprintf(f,'NODE%dX\t',i);
   if cnd(clnk(i,1))
      fprintf(f,'NODE%dX\t',cnd(clnk(i,1)));
   else
      fprintf(f,'ARRY%dX\t',cidx(clnk(i,1)));
   end
   if cnd(clnk(i,3))
      fprintf(f,'NODE%dX\t',cnd(clnk(i,3)));
   else
      fprintf(f,'ARRY%dX\t',cidx(clnk(i,3)));
   end
   cnd(clnk(i,1))=i;
   fprintf(f,'%f%s',(clnk(i,5)-min(clnk(:,5)))/(max(clnk(:,5)+1)-min(clnk(:,5)))*1.8-0.8,char([13 10]));
end
fclose(f);


%ls=cell((N+2),1);
%f2=fopen(ofname,'r');
%for i=1:(N+2)
%   l=fgetl(f2);
%   ls{i}=splitstr(l,[9]);
%end
%fclose(f2);

ls=read_dlm_file(ofname);
% is it a GCT?
if strcmp(ls{1}{1},'#1.2') 
  ls=ls(3:end);
end

% assumes no EWEIGHT or GWEIGHT
f=fopen([fname '.cdt'],'w');
fprintf(f,'%s\t%s\t%s\t%s\t','GID',ls{1}{1},ls{1}{2},'GWEIGHT');
for j=1:(cN-1)
   fprintf(f,'%s\t',ls{1}{cidx(j)+2});
end
fprintf(f,'%s%s',ls{1}{cidx(cN)+2},char([13 10]));

fprintf(f,'%s\t\t\t\t','AID');
for j=1:(cN-1)
   fprintf(f,'ARRY%dX\t',cidx(j));
end
fprintf(f,'ARRY%dX%s',cidx(end),char([13 10]));

fprintf(f,'%s\t','EWEIGHT');
fprintf(f,'\t\t\t');
for j=1:(cN-1)
   fprintf(f,'%s\t','1');
end
fprintf(f,'%s%s','1',char([13 10]));

if strcmp(ls{2}{1},'EWEIGHT')
  skip_2nd=1;
else
  skip_2nd=0;
end

for i=1:N
   fprintf(f,'GENE%dX\t%s\t%s\t%s\t',...
	   idx(i),...
	   ls{idx(i)+1+skip_2nd}{1},...
	   ls{idx(i)+1+skip_2nd}{2},...
	   '1');
	for j=1:(cN-1)
	  fprintf(f,'%s\t',ls{idx(i)+1+skip_2nd}{cidx(j)+2});
	end
	fprintf(f,'%s%s',ls{idx(i)+1+skip_2nd}{cidx(cN)+2},char([13 10]));
end
fclose(f);
