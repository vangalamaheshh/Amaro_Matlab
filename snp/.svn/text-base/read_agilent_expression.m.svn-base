function D=read_agilent_expression(fname)

[f,fid]=read_dlm_file(fname,char(9),1);
form=[ '%s%s%s%s' repmat('%f',1,length(f{1})-4)];
F=textscan(fid,form,'delimiter',char(9),'bufSize',1000000);

D.sdesc=f{1}(5:end);
D.dat=cat(2,F{5:end});
D.gacc=F{3};
D.gdesc=F{4};
D.id=F{2};
D.symb=F{1};
fclose(fid);
