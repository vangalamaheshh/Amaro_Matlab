function transpose_table(in_fname,out_fname)

x=read_dlm_file(in_fname);
x=cat(1,x{:});


st=[repmat('%s\t',1,size(x,1)-1) '%s\n'];

f=fopen(out_fname,'w');
fprintf(f,st,x{:});
fclose(f);

