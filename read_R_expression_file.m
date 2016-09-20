function D=read_R_expression_file(fname)

f=fopen(fname,'r');
l=read_dlm_file(f,char(9),1);
l=l{1};
nsamples=length(l);
dat=textscan(f,['%s' repmat('%f',1,nsamples)]);
fclose(f);

D.sdesc=l;
D.gacc=dat{1};
D.gdesc=D.gacc;
D.dat=cat(2,dat{2:end});
