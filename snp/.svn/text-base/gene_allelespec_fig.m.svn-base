function gene_allelespec_fig(genelist,infile, mutfile, ext,n,refgene)
%
%gene_allelespec_fig(genelist,infile, mutfile, outfile,n, refgene)
%
%genelist=text file of genes, one gene per row
%infile=struct containing min and max info in .adat fields(.mat file)-needs either linkingid in .sis or TCGA names,
%ext=unique name for output table of mutation info and gene plots, 
%n=#of kb window on either side of gene cds region, if there are no SNPs in cds region
%refgene .mat file, if not hg18 - use [] if hg18

if ~exist(refgene)
refgene_file  ='/xchip/tcga/Annotation_Files/UCSC_hg18/hg18_20070817.mat';
load(refgene_file);
else
    load(refgene)
end
rg=add_chrn(rg);
rg=mark_futreal_genes(rg);
for i=1:length(rg)
    rg(i).symbol=rg(i).symb;
end
rg1=collapse_rg_to_genes(rg);
rgs=regexprep({rg1.symbol},'[\*]+','');

if ~exist(infile)
    load('Affy128newseg.mat')
else
    x=load(infile)
    fn=fieldnames(x);
    D=getfield(x,fn{1})
    Dmin=D;
    Dmin.dat=D.adat(:,:,1);
    Dmax=Dmin;
    Dmax.dat=D.adat(:,:,2);
end

list=read_dlm_file(genelist);
d=cat(1,list{:});
[Mt,m1,m2]=match_string_sets_hash(d(:,1),rgs);
rgg=rg1(m2);

mut=read_dlm_file(mutfile);
% mut=read_dlm_file('/xchip/tcga/gbm/analysis/Barbara/TCGA_GBM_Level3_All_Somatic_Mutations_DataFreeze2.maf') % TCGA list
% mut=read_dlm_file('/xchip/cancergenome/cancergenome02/data02/Barbara/mutation_files/TSP_20080515_hg17.mut') % TSP list
K=[];
clear k j
f=regexp(mutfile,'.maf','match')
if strmatch(f,'.mut') %%%% for .mut file
    for i=1:length(mut)
        k(i)=mut{i}(4);
        j(i)=mut{i}(6);
    end
    K={k j};
    
elseif strmatch(f,'.maf') %%% for .maf file
    for i=1:length(mut)
        k(i)=mut{i}(16);
        j(i)=mut{i}(1);
        s(i)=mut{i}(9);
    end
    
    k1=regexp(k,'T(\w*)-(\w*)-(\w*)', 'match');
    k2={''};
    for i=1:length(mut)
        k2=[k2 k1{i}];
    end
    
    K2={k2 j s};
    ids=strmatch('Silent_Mutation',K2{3});
    K{1}=K2{1}(setdiff(1:length(K2{3}),ids));
    K{2}=K2{2}(setdiff(1:length(K2{3}),ids));
    K{3}=K2{3}(setdiff(1:length(K2{3}),ids));
end

if isfield(Dmin.sis, 'linkingid')
    D_info={Dmin.sis.linkingid};  %tsp
else
    D_info=regexp({Dmin.sis.name},'T(\w*)-(\w*)-(\w*)', 'match');  %tcga only
end
D_info=cat(1,D_info{:});
[KD,idk,idd]=match_string_sets_hash(K{1},D_info);
K1{1}=K{1}(idk);
K1{2}=K{2}(idk);
K1{3}={Dmin.sis(idd).name};
J=struct('name',K1{3},'mut',K1{2});
X=struct('name',unique({J.name}));

for i=1:length(X)
    X(i).mut={J(strmatch(X(i).name,{J.name})).mut};
end
[XJ,idx,idd]=match_string_sets({X.name},{Dmin.sis.name});
mutations={X(idx).mut};
Dmin.sdesc={Dmin.sis.array};
for i=1:length(idd)
    Dmin.sis(idd(i)).mut=mutations(i);
end
name=[ext '_mutation_summary.txt'];
f=fopen(name,'w');
fprintf(f,'%s\t%s\t%s\n','gene','# mutated','# sequenced')
for i=1:length(rgg)
    %subplot((length(d)*10),1,(((10*i)-5):10*i))
    figure;
    h3=[];
    if nansum(nanmean(Dmin.dat(Dmin.chrn==rgg(i).chrn & Dmin.pos>rgg(i).cds_start & Dmin.pos<rgg(i).cds_end,:))) ~=0
        m=0;
    else
        m=n;
    end
    dmin=nanmean(Dmin.dat(Dmin.chrn==rgg(i).chrn & Dmin.pos>rgg(i).cds_start-(m*1000) & Dmin.pos<rgg(i).cds_end+(m*1000),:));
    dmax=nanmean(Dmax.dat(Dmax.chrn==rgg(i).chrn & Dmax.pos>rgg(i).cds_start-(m*1000) & Dmax.pos<rgg(i).cds_end+(m*1000),:));
    if nansum(dmin)~=0
        h1=scatter(dmin,dmax,8,'gx');
        hold on
        c=0; t=0;
        for j=1:length({Dmin.sis.mut})
            if ~isempty(Dmin.sis(j).mut)
                idm=strmatch({rgg(i).symbol},Dmin.sis(j).mut{1});
                mutmin=nanmean(Dmin.dat(Dmin.chrn==rgg(i).chrn & Dmin.pos>rgg(i).cds_start-(m*1000) & Dmin.pos<rgg(i).cds_end+(m*1000),j));
                mutmax=nanmean(Dmax.dat(Dmax.chrn==rgg(i).chrn & Dmax.pos>rgg(i).cds_start-(m*1000) & Dmax.pos<rgg(i).cds_end+(m*1000),j));
                t=t+1;
                if isempty(idm)
                    h2=scatter(mutmin,mutmax,20,'bx');
                else
                    h3=scatter(mutmin,mutmax,20,'r', 'filled');
                    c=c+1;
                end
            end
        end
        fprintf(f,'%s\t%d\t%d\n',rgg(i).symbol,c,t)
        display(rgg(i).symbol);
        display([c t]);
        add_x_equ_y_line;
        title(rgg(i).symb)
        ylabel('max')
        xlabel('min')
        hAnnotation = get(h1,'Annotation');
        set(get(hAnnotation','LegendInformation'),'IconDisplayStyle','on')
        set(h1,{'DisplayName'},{'not sequenced'}')
        set(h2,{'DisplayName'},{'sequenced, not mutated'}')
%         if ~isempty(h3)
            set(h3,{'DisplayName'},{'mutated'}')
%         else
%             h3=[];
%         end
        legend([h1 h2 h3],4)
        print_D([ext '_MinMax_' rgg(i).symbol],{{'pdf'},{'png','-r180'}});
    end
end
fclose(f);




