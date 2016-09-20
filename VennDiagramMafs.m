

MAFS=load_struct('/xchip/cga_home/amaro/Brain_Mets/maf_comparison/maf_mapping.txt');

for i=1:22
v{i}=num2str(i);
end
ffg={};
for i=1:slength(MAFS)
maf_new=load_struct(MAFS.new_maf{i});
maf_new=trimStruct(maf_new,ismember(maf_new.Chromosome,v));

    cols=[1:10 80:82];
    call_stats=load_table(MAFS.call_stats{i},[],[],cols);%cols);
    
    

maf_old=load_struct(MAFS.old_maf{i});

for j=1:slength(maf_new)
maf_new.site{j,1}=sprintf('%s-%s',maf_new.Chromosome{j},maf_new.Start_position{j});
end

for j=1:slength(maf_old)
maf_old.site{j,1}=sprintf('%s-%s',maf_old.Chromosome{j},maf_old.Start_position{j});
end

for j=1:length(call_stats.contig)
call_stats.site{j,1}=sprintf('%s-%d',call_stats.contig{j},call_stats.position(j));
end

venn_counts.Both(i,1)=sum(ismember(maf_old.site,maf_new.site));
venn_counts.Old(i,1)=sum(~ismember(maf_old.site,maf_new.site));
venn_counts.New(i,1)=sum(~ismember(maf_new.site,maf_old.site));
venn_counts.pair_id{i,1}=MAFS.pair_id{i,1};




sti=~ismember(maf_old.site,maf_new.site);
sites={maf_old.site{sti}};

NumberOfSitesMissed(i)=venn_counts.Old(i)-sum(ismember(call_stats.site,sites));




failure_reasons={call_stats.failure_reasons{ismember(call_stats.site,sites)}}';
ffg=vertcat(ffg,failure_reasons);


end