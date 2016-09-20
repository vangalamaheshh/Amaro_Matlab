function run_R_file(R_call_fname,R_code_fname,R_func,params_str)

f=fopen(R_call_fname,'w');
fprintf(f,['source("' R_code_fname '")' newline]);

fprintf(f,[R_func '(' params_str ')' newline ]);
fprintf(f,['q()' newline]);
fclose(f);

