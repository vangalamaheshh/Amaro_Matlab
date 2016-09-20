% cd /xchip/cga_home/amaro/TumorInNormal/Pancan_NT_paper/PanCanPsetDance/
% while read p; do
%   pset=$(fiss pset_list $p | head -1)
%   fiss annot_get $p pair pset=$pset call_stats_capture > $p.call_stats.tsv
%   fiss annot_get $p sample pset=$pset clean_bam_file_capture short_letter_code sample_type recapseg_seg_file > $p.tsv
% done <Pran_TCGA_List



STSVs=load_struct_noheader('~/Projects/TumorInNormal/Pancan_NT_paper/ListOfSampleTables.tsv');
C=struct;
for s=1:slength(STSVs)
    
    S=load_struct(STSVs.col1{s});
    STSVs.col1{s}
    for i=1:slength(S)
        if s>13
            S.indiv{i,1}=S.sample_id{i}(1:17);
            S.Ttype{i,1}=S.sample_id{i}(1:4);
        else
            S.indiv{i,1}=S.sample_id{i}(1:12);
            strs=regexp(STSVs.col1{s},'/','split');
            str=regexp(strs{end},'_|.tsv','split');
            Ttype=str{end-1};
            S.sample_id{i}=strcat([Ttype,'_',S.sample_id{i}]);
            S.Ttype{i,1}=Ttype;
        end
    end
    indivs=unique(S.indiv);
    for j=1:length(indivs)
        if sum(ismember(S.short_letter_code(ismember(S.indiv,indivs{j})),{'NT'}))>=1
            S.KEEP(ismember(S.indiv,indivs{j}),1)={'YES'};
        else
            S.KEEP(ismember(S.indiv,indivs{j}),1)={'NO'};
        end
    end
    S_keep=reorder_struct(S,ismember(S.KEEP,'YES'));
    C=mergeStruct(S_keep,C);
    
    
end
    



PTSVs=load_struct_noheader('~/Projects/TumorInNormal/Pancan_NT_paper/ListOfPairTables.tsv');
PC=struct;
for pp=1:slength(PTSVs)
    P=load_struct(PTSVs.col1{pp});
    PTSVs.col1{pp}
    P.KEEP=repmat({'NO'},slength(P),1);
    for p=1:slength(P)
        for i=1:slength(C)
            if ~isempty(strfind(P.pair_id{p},C.indiv{i}))
                P.KEEP{p,1}='YES';
                P.Ttype{p,1}=C.Ttype{i};
   
                break
            end
        end
    end
    P_keep=reorder_struct(P,ismember(P.KEEP,'YES'));
    if isfield(P_keep,'Ttype')
    PC=mergeStruct(P_keep,PC);
    count(PC.Ttype)
    end
   
end
PC=rmfield(PC,'N');
PC=reorder_struct(PC,cellfun(@isempty,strfind(PC.pair_id,'NB-NT')));
PC=reorder_struct(PC,cellfun(@isempty,strfind(PC.pair_id,'NT-NB')));
PC=reorder_struct(PC,cellfun(@isempty,strfind(PC.pair_id,'TP-TP')));
PC=reorder_struct(PC,cellfun(@isempty,strfind(PC.pair_id,'NT-TP')));
PC=reorder_struct(PC,cellfun(@isempty,strfind(PC.pair_id,'NB-TP')));

for i=1:slength(PC)
strs=regexp(PC.pair_id{i},'NT|NB|NWGA','match');
if ~isequal(strs,{})
PC.NormalType{i,1}=strs{1};
else    
PC.NormalType{i,1}='NaN';  
end
end
PC=reorder_struct(PC,~ismember(PC.NormalType,'NaN'));
save_struct(PC,'~/Projects/TumorInNormal/Pancan_NT_paper/PanCanPairsForUpload.tsv');
save_struct(C,'~/Projects/TumorInNormal/Pancan_NT_paper/PanCanSamplesForUpload.tsv');


























