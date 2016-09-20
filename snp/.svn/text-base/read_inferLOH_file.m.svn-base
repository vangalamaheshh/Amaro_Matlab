function LOH=read_inferLOH_file(fname)

if (1)
unix([ 'head -2 ' fname ' > header.txt' ]);
unix([ 'tail +3 ' fname ' | sed ''s/L\t/1\t/g''  > tmp1.txt']);
unix([ 'sed ''s/R\t/2\t/g'' tmp1.txt > tmp2.txt']);
unix([ 'sed ''s/N\t/3\t/g'' tmp2.txt > tmp3.txt']);
unix([ 'perl -e "@k=<>;splice @k, -6 if @k > 6; print @k" < tmp3.txt > tmp.txt']);
unix([ 'cat header.txt tmp.txt > tmp4.txt']);
unix([ 'rm header.txt tmp.txt tmp1.txt tmp2.txt tmp3.txt']);
end

% L==1, R==2, N==3

[f,fid]=read_dlm_file('tmp4.txt',char(9),2);
form=[ '%s%s%s%f%f%f'  repmat('%f',1,length(f{2})-6)];
LOH_dat=textscan(fid,form);
LOH.marker=LOH_dat{1};
%LOH.dbSNP=LOH_dat{2};
LOH.chr=LOH_dat{3};
LOH.pos=LOH_dat{4}*1.0e06;
LOH.cM=LOH_dat{5};
LOH.score=LOH_dat{6};
LOH.dat=cat(2,LOH_dat{7:end-1}); 

LOH.sdesc=f{2}(7:end-1);

%f2=read_dlm_file(fid,char(9),1);

%if ~isempty(f2)
%  keyboard
%  LOH=reorder_D_rows(LOH,1:(size(LOH.dat,1)-3));
%  LOH=rmfield(LOH,{'orig','origidx','gorigidx','history'});
%end

fclose(fid);
