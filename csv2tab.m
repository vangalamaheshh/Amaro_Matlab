function csv2tab(infname,outfname)

fin=fopen(infname,'r');
fout=fopen(outfname,'w');
lnum=0;
while(1)
  tline = fgetl(fin);
  if ~ischar(tline), break, end
  qpos=find(tline=='"');
  if isempty(qpos)
      ln=dlmsep(tline,',');
      fprintf(fout,['%s' repmat('\t%s',1,length(ln)-1) '\n'],ln{:});
  else
      tline=regexprep(tline,'""',char([1 0])); %[remove, turn to "]
      qpos=find(tline=='"');
      if mod(length(qpos),2)~=0
          error('no equal number of "');
      end
      for i=1:2:length(qpos)
          cpos=find(tline(qpos(i):qpos(i+1))==',');
          tline(qpos(i)+cpos-1)=2;
          tline(qpos([i i+1]))=char(1); % remove
      end   
      tline(tline==char(0))='"';
      tline(tline==char(1))=[];
      ln=dlmsep(tline,',');
      for j=1:length(ln)
          ln{j}(ln{j}==char(2))=',';
      end
       fprintf(fout,['%s' repmat('\t%s',1,length(ln)-1) '\n'],ln{:});     
  end   
  lnum=lnum+1;
  if mod(lnum,1000)==0
    disp(lnum)
  end
end
fclose(fin);
fclose(fout);
