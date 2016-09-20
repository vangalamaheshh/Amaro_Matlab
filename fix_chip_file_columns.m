function fix_chip_file_columns(infname,outfname,chip_name)
c={'Probe','Symbol','Title'};

fin=fopen(infname,'r');
ln=fgetl(fin);
tabs=dlmsep(ln,char(9));
tabs(find(cat(1,cellfun('isempty',tabs))))=[];
cl_st=[];
cl=[];
for j=1:length(c)
    tmp=find(~cat(1,cellfun('isempty',regexpi(tabs,lower(c{j})))));
    if ~isempty(tmp)
        cl_st=[cl_st ',$' num2str(tmp(1))];
        cl(j)=tmp(1);
    else
        cl_st=[cl_st ',"NA"'];
        cl(j)=NaN;
    end
end

if length(tabs)>length(c)
    cl_st=[cl_st sprintf(',$%d',setdiff(1:length(tabs),cl))];
end 

unix_st=[ 'tail +2 ' infname ' | awk -F''\t'' ''{ printf "' repmat('%s\t',1,length(tabs)) '%s\n",' ...
        '"' chip_name '"' cl_st '}'' > ' outfname ];
disp(unix_st);
unix(unix_st);
