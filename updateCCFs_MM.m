MM.key=strcat(MM.sample,MM.Start_position,MM.Tumor_Seq_Allele2);
mafs=load_struct('~/Projects/CLL_Rush/ABSOLUTE/ABS_mafs.txt');

for i=1:slength(mafs)
m=load_struct(mafs.mafs{i});
fn=fieldnames(m);
m.key=strcat(m.sample,m.Start_position,m.Tumor_Seq_Allele2);

    for j=1:slength(m)
     l=find(ismember(MM.key,m.key{j}),1);
    if ~isempty(l)
     for f=1:length(fn)
     MM.(fn{f}){l}=m.(fn{f}){j};
     end
    end

    end
end