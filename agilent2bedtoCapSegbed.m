fid=fopen('~/Projects/Cranios/WESv2/S0293689_Covered.bed')
i=1;

while ischar(tline)
     
long=regexp(tline,'\t','split');
chr=regexp(long{1},'chr','split');                
chromosome{i}=chr{2};                          
start{i}=long{2};
stop{i}=long{3};
tar_str=regexp(long{4},',','split');           
for j=1:size(tar_str,2)
if ~isempty(regexp(tar_str{j},'ref','match')) 
ss=regexp(tar_str{j},'ref|','split')
gene=ss{2};
break
end
end
target{i}=sprintf('%s-bait%s',gene,num2str(i));
i=i+1;
tline = fgetl(fid);
end
 
fid2=fopen('~/Projects/Cranios/WESv2/S0293689_Covered.csv','w')l
for i=1:size(chromosome,2)
fprintf(

end