x=load_struct('/xchip/cga_home/amaro/CLL/ABSOLUTE/TimePoint1.mafs');

M=load_struct(x.maf{1});

for i=2:slength(x)
m=load_struct(x.maf{i});
M=mergeStruct(M,m);
M=rmfield(M,'N');
if mod(i,5)==0
    i
end
end
    

save_struct(M,'/xchip/cga_home/amaro/CLL/ABSOLUTE/CombinedABSResultsForDanTP1.maf')
clear M m x

x=load_struct('/xchip/cga_home/amaro/CLL/ABSOLUTE/TimePoint2.mafs');

M=load_struct(x.maf{1});

for i=2:slength(x)
m=load_struct(x.maf{i});
M=mergeStruct(M,m);
M=rmfield(M,'N');
if mod(i,5)==0
    i
end
end
    

save_struct(M,'/xchip/cga_home/amaro/CLL/ABSOLUTE/CombinedABSResultsForDanTPLater.maf')
