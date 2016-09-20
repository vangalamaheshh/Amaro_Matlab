function test_gpm(fname,inp2)

f=fopen(fname,'w');

fprintf(f,'%f\n',str2num(inp2));
fclose(f);

% quit - used in the past

