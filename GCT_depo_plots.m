data=load_table('~/Documents/GCT_mitochondria_data/AggregateGCTdata.txt');
tumors={'T1','T2','T3','T4','T5','T6','T7'};
normals={'N1','N1','N2','N3','N4','N5','N6'}; %N1 has two tumors (T1,T2)
peptides=fieldnames(data);
peptides={peptides{4:end}};

peptides={'Bim_100','Bad_100','Hrk_100','Puma_100'};
figure()
hold on

for i=1:length(peptides)
    for j=1:length(tumors)
        kmt=find(ismember(data.sample_type,strcat([tumors{j},'_mean'])));
        ket=find(ismember(data.sample_type,strcat([tumors{j},'_err'])));
        kmn=find(ismember(data.sample_type,strcat([normals{j},'_mean'])));
        ken=find(ismember(data.sample_type,strcat([normals{j},'_err'])));
        
        errorbar([i-.25,i+.25],[data.(peptides{i})(kmn),data.(peptides{i})(kmt)],[data.(peptides{i})(ken),data.(peptides{i})(ket)],'r','Marker','.','MarkerSize',1,'LineWidth',2,'LineStyle','-')
        
        if data.(peptides{i})(kmn)>data.(peptides{i})(kmt)
            disp(sprintf('This tumor bucks the trend : %s with peptide %s',tumors{j},peptides{i}))
        end
    end
end

for i=1:length(peptides)
    for j=1:length(tumors)
        kmt=find(ismember(data.sample_type,strcat([tumors{j},'_mean'])));
        ket=find(ismember(data.sample_type,strcat([tumors{j},'_err'])));
        kmn=find(ismember(data.sample_type,strcat([normals{j},'_mean'])));
        ken=find(ismember(data.sample_type,strcat([normals{j},'_err'])));
        
        normals_dat(i,j)=data.(peptides{i})(kmn);
        tumors_dat(i,j)=data.(peptides{i})(kmt);
        plot([i-.25],[data.(peptides{i})(kmn)],'.','MarkerSize',30,'Color',[.5 .5 .5])
        plot([i+.25],[data.(peptides{i})(kmt)],'.','MarkerSize',30,'Color',[0 0 0])
        
    end
end

peptides={'Bim','Bad','Hrk','Puma'};
ylim([0,1])
ylabel('% Depolarization','FontSize',30)
set(gca,'XTick',[1,2,3,4],'XTickLabel',peptides,'FontSize',20)



peptides={'DMSO','Bim_100','Bad_100','Hrk_100','Puma_100'};



