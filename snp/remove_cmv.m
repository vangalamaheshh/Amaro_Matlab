function C=remove_cmv(C,cmv_file)

cmv=tab=read_table(cmv_file,'%s%s',char(9),0);
read_table('%d\t%d\n',
