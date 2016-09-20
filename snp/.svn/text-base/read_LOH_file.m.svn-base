function LOH=read_LOH_file(fname)

if (1)
unix([ 'head -2 ' fname ' > header.txt' ]);
unix([ 'tail +3 ' fname ' | sed ''s/L\t/1\t/g''  > tmp1.txt']);
unix([ 'sed ''s/R\t/2\t/g'' tmp1.txt > tmp2.txt']);
unix([ 'sed ''s/N\t/3\t/g'' tmp2.txt > tmp3.txt']);
unix([ 'cat header.txt tmp3.txt > tmp4.txt']);
end

% L==1, R==2, N==3

[f,fid]=read_dlm_file('tmp4.txt',char(9),2);
form=[ '%s%s%f%f%s'  repmat('%f',1,length(f{2})-5)];
LOH_dat=textscan(fid,form);
LOH.marker=LOH_dat{1};
LOH.chr=LOH_dat{2};
LOH.pos=LOH_dat{3};
LOH.cM=LOH_dat{4};
LOH.score=LOH_dat{5};
LOH.dat=cat(2,LOH_dat{6:(end-1)}); 

LOH.sdesc=f{2}(6:(end-1));

f2=read_dlm_file(fid,char(9),1);

if ~isempty(f2)
  keyboard
  LOH=reorder_D_rows(LOH,1:(size(LOH.dat,1)-3));
  LOH=rmfield(LOH,{'orig','origidx','gorigidx','history'});
end

fclose(fid);
       
 
