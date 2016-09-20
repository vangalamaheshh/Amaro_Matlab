function res=file_exist(fname)

f=fopen(fname,'r');
if f < 0 
    res=0;
else
    fclose(f);
    res=1;
end
