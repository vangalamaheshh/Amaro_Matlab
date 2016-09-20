function [ch,st,en]=str2genomic_location(str)

col=find(str==':');
ch=chromosome2num(str(1:(col(1)-1)));
minus=find(str=='-');
st=str2num(str((col(1)+1):(minus(1)-1)));
en=str2num(str((minus(1)+1):end));

