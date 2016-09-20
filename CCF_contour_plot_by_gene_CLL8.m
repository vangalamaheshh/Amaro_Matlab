sigABS=load_table('/Users/amaro/Documents/StilgenbauerFigureMafforCCFprogressionbyGene.maf');
for i=1:slength(sigABS)
    sigABS.indvidual{i,1}=sigABS.sample{i}(1:13);
end

sigABS.key=strcat(num2str(sigABS.Start_position),sigABS.indvidual);

genes=unique(sigABS.Hugo_Symbol);
for s=1:1%length(genes)
    genes{s}='TP53';
    all=[];
   g_maf=reorder_struct(sigABS,ismember(sigABS.Hugo_Symbol,genes{s})); 
   hold on   


[n l]=count(sigABS.key(ismember(sigABS.Hugo_Symbol,genes{s})&ismember(sigABS.TP,1)));
    connect=l(n>1);
    all.key=l;
    all.map=zeros(slength(all),1);
    map=zeros(101,101);
        f_g_maf=fieldnames(g_maf);
    ccff_start=find(ismember(f_g_maf,'x0'));
    ccf_fields=f_g_maf(ccff_start:ccff_start+100);
for i=1:slength(all)
    

    TPs=unique(g_maf.TP(ismember(g_maf.key,all.key{i})));
    if ismember(2,TPs)
    all.ccf(i,1)=g_maf.ccf_hat(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    all.ccf_high(i,1)=g_maf.ccf_CI95_high(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    all.ccf_low(i,1)=g_maf.ccf_CI95_low(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1));
    range1=[all.ccf_low(i):.01:all.ccf_high(i)]';

    
    
    all.ccf2(i,1)=g_maf.ccf_hat(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));
    all.ccf2_high(i,1)=g_maf.ccf_CI95_high(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));
    all.ccf2_low(i,1)=g_maf.ccf_CI95_low(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2));
    range2=[all.ccf2_low(i):.01:all.ccf2_high(i)]';
    if all.ccf2(i)>0||all.ccf(i)>0
    all.map(i,1)=1;    
    for j=1:length(range1)
        
        index_j=round(range1(j)*100)+1;
        ccf_field_j=ccf_fields{index_j};
        for k=1:length(range2)
            
            index_k=round(range2(k)*100)+1;
            ccf_field_k=ccf_fields{index_k};
            map(index_k,index_j)=g_maf.(ccf_field_k)(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,2))+g_maf.(ccf_field_j)(ismember(g_maf.Hugo_Symbol,genes{s})&ismember(g_maf.key,all.key{i})&ismember(g_maf.TP,1))+map(index_j,index_k);

        end
    end
    end
    
    else
        all.ccf2(i,1)=NaN;
    end
end

all=reorder_struct(all,~isnan(all.ccf2));


end