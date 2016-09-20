maf='~/Documents/SignificantGeneMaf.CLL8.tsv';
sig_gene_list='/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/1_09_2015_PoNCut4_WithSalvage/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt';
M=load_struct(maf); % A Maf with Standard columns
S=load_struct(sig_gene_list); % has to have gene and q columns

uni_prot_mapping=load_struct('/Users/amaro/Downloads/UniProtLengthGeneIDmap'); %This is the result of a query to UNIPROT: http://www.uniprot.org/
cosmic_track=load_struct('/Users/amaro/Downloads/CosmicTrackStickFigs.tsv');




S.q=str2double(S.q);
S=reorder_struct(S,S.q<.1);


for i=1:slength(uni_prot_mapping)
    str=split(uni_prot_mapping.Gene_names{i},' ');
    uni_prot_mapping.Hugo{i,1}=str{1};
end
M=reorder_struct(M,ismember(M.Hugo_Symbol,S.gene)); 
M=reorder_struct(M,ismember(M.Hugo_Symbol,'RPS15'));
%M.Protein_Change(ismember(M.Hugo_Symbol,'MYD88')&ismember(M.Protein_Change,'p.*160R'))={'p.L265P'};
cosmic_track=reorder_struct(cosmic_track,ismember(cosmic_track.Genename,S.gene));
%Headers required to go into mutfig
%Hugo_Symbol	Feature	label	Chromosome	type	start	end	pair_id	color	Grouping	groupchar	gene
S=reorder_struct(S,ismember(S.gene,'RPS15'));
cosmic_trans='ENST00000233609';
for i=1:slength(M)
    if ~isequal(M.Other_Transcripts{i},'')
    transcripts=split(M.Other_Transcripts{i},'|');
    a=strfind(transcripts,'ENST00000233609');
    new_trans=transcripts(~cellfun(@isempty,a));
    new_protien_change=regexp(new_trans,'p.\S+','match');
    M.Protein_Change(i)=new_protien_change{1};
    end
end

for i=1:slength(S)
m.Hugo_Symbol=M.Hugo_Symbol(ismember(M.Hugo_Symbol,S.gene{i}));
m.feature=repmat({'mutation'},slength(m),1);
m.label=M.Protein_Change(ismember(M.Hugo_Symbol,S.gene{i}));
m.Chromosome=M.Chromosome(ismember(M.Hugo_Symbol,S.gene{i}));
m.type=tr(M.Variant_Classification(ismember(M.Hugo_Symbol,S.gene{i})),'_',' ');
m.start=M.UniProt_AApos(ismember(M.Hugo_Symbol,S.gene{i}));
m.end=M.UniProt_AApos(ismember(M.Hugo_Symbol,S.gene{i}));
m.pair_id=M.sample(ismember(M.Hugo_Symbol,S.gene{i}));
m.Grouping=repmat({'3'},slength(m),1);
m.groupchar=m.type;
m.gene=m.Hugo_Symbol;
m.Hugo_Symbol{end+1}=m.Hugo_Symbol{1};
m.feature{end+1}='protein';
m.label{end+1}=m.Hugo_Symbol{1};
m.Chromosome{end+1}=m.Chromosome{1};
m.type{end+1}='NA';
m.start{end+1}='1';
m.end(end+1)=uni_prot_mapping.Length(ismember(uni_prot_mapping.Hugo,S.gene{i}));
m.pair_id{end+1}='NA';
m.Grouping{end+1}='1';
m.groupchar{end+1}='Other';
m.gene{end+1}=m.Hugo_Symbol{1};
gene_cosmic=reorder_struct(cosmic_track,ismember(cosmic_track.Genename,S.gene{i}));
[mn nm]=count(gene_cosmic.MutationAA);
sites_to_add=nm(mn>10);

for j=1:length(sites_to_add)
strs=regexp(sites_to_add{j},'p.\S','split');
if ~isempty(strs{2})
strs=regexp(strs{2},'[A-Z]?','split');
end



m.start{end+1}=strs{1};
m.end{end+1}=strs{1};
m.type{end+1}='Cosmic';
m.pair_id{end+1}='NA';
m.feature{end+1}='Domain';
m.Grouping{end+1}='1';
m.groupchar{end+1}='Other';
m.gene{end+1}=m.gene{1};
m.label{end+1}=strs{1};
m.Hugo_Symbol{end+1}=m.Hugo_Symbol{1};
m.Chromosome{end+1}=m.Chromosome{1};
end
for ii=1:slength(m)
strs=regexp(m.label{ii},'p.[*A-Za-z-_]+','split');
if length(strs)>1 && ~isequal(m.Hugo_Symbol{ii},m.label{ii}) && isequal(m.feature{ii},'mutation')
strs=regexp(strs{2},'[?A-Z*?a-z_?]+','split');
m.start{ii,1}=strs{1};
m.end{ii,1}=strs{1};
end
end


m.color=zeros(slength(m),3);


m.color(ismember(m.type,'3''UTR'),:)=repmat([166/255,206/255,227/255],sum(ismember(m.type,'3''UTR')),1);
m.color(ismember(m.type,'5''UTR'),:)=repmat([166/255,206/255,227/255],sum(ismember(m.type,'5''UTR')),1);
m.color(ismember(m.type,'De novo Start OutOfFrame'),:)=repmat([.5,.5,.5],sum(ismember(m.type,'De novo Start OutOfFrame')),1);
m.color(ismember(m.type,'Frame Shift Del')|ismember(m.type,'Frame Shift Ins'),:)=repmat([255/255,127/255,0/255],sum(ismember(m.type,'Frame Shift Del')|ismember(m.type,'Frame Shift Ins')),1);
m.color(ismember(m.type,'In Frame Del')|ismember(m.type,'In Frame Ins'),:)=repmat([253/255,191/255,111/255],sum(ismember(m.type,'In Frame Del')|ismember(m.type,'In Frame Ins')),1);
m.color(ismember(m.type,'Nonsense Mutation')|ismember(m.type,'Nonstop Mutation'),:)=repmat([227/255,26/255,28/255],sum(ismember(m.type,'Nonsense Mutation')|ismember(m.type,'Nonstop Mutation')),1);
m.color(ismember(m.type,'Missense Mutation'),:)=repmat([95/255,95/255,44/255],sum(ismember(m.type,'Missense Mutation')),1);
m.color(ismember(m.type,'Splice Site'),:)=repmat([202/255,178/255,214/255],sum(ismember(m.type,'Splice Site')),1);
m.color((ismember(m.type,'Cosmic')),1:3)=repmat([1 0 0],sum(ismember(m.type,'Cosmic')),1);
m.color((ismember(m.type,'NA')),1:3)=repmat([.75 .75 .75],sum(ismember(m.type,'NA')),1);
m.color((ismember(m.type,'Silent')),1:3)=repmat([.75 .75 .75],sum(ismember(m.type,'Silent')),1);
m.color((ismember(m.type,'Intron')),1:3)=repmat([.75 .75 .75],sum(ismember(m.type,'Intron')),1);
m.color((ismember(m.type,'Start Codon SNP')),1:3)=repmat([.75 .75 .75],sum(ismember(m.type,'Start Codon SNP')),1);

m
m=reorder_struct(m,~ismember(m.label,''));

mutfig(m,strcat('~/Documents/RevisionFiguresCLL8/StickFigure/',S.gene{i},'_.png'))
close all
clear m
end



