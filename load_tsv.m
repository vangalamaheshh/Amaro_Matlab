function X = load_tsv(tsv_file,comment,delimiter)
% load_tsv(tsv_file,OPT) skips header lines 
if nargin<2
    comment='#';
end
if nargin<3
    delimiter='\t';
end
cmd =['grep "^' comment '" ' tsv_file ' | wc -l'];
[o nh]=unix(cmd);
nh=round(abs(str2num(nh)));
if isempty(nh), nh=0;end
X0 = readtable(tsv_file,'Delimiter',delimiter,'FileType','text','HeaderLines',nh);
X = table2struct(X0,'ToScalar',true);

