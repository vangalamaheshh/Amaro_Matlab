function name=index2name(idx)

if (size(idx,1)==1)
    idx=idx';
end
name=cell(size(idx));
name(idx>22)=cellstr(char(idx(idx>22)+'A'));
name(idx<=22)=cellstr(num2str(idx(idx<=22)));
