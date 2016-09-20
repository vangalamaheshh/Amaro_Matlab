function [ amplicons ] = MAFtoAmplicons( Maf,sig_genes,type )
%Takes in a Maf and significant genes file and produces a list of amplicons
%covering the mutations in the MAF
S=load_struct(sig_genes);
S.q=cellfun(@str2num,S.q);
S=reorder_struct(S,S.q<.1);


if isequal(type,'STAD')
load('~/Projects/STAD/RPKM_stats.mat');
for i=1:slength(S)
    if isempty(RPKM_stats.median(ismember(RPKM_stats.gene,S.gene{i})))
    S.RPKM=0;
    else
    S.RPKM(i,1)=cell2mat(RPKM_stats.median(find(ismember(RPKM_stats.gene,S.gene{i})==1,1)));
    end
end
S=reorder_struct(S,S.RPKM>1);
end
M=load_struct(Maf); 
M=reorder_struct(M,ismember(M.Hugo_Symbol,S.gene));
M=reorder_struct(M,ismember(M.Variant_Classification,'In_Frame_Del')|ismember(M.Variant_Classification,'In_Frame_Ins')|...
    ismember(M.Variant_Classification,'Missense_Mutation')|ismember(M.Variant_Classification,'Nonsense_Mutation')|...
    ismember(M.Variant_Classification,'Splice_Site')|ismember(M.Variant_Classification,'Stop_Codon_Ins')|...
    ismember(M.Variant_Classification,'Translation_Start_Site')|ismember(M.Variant_Classification,'Frame_Shift_Del')|...
    ismember(M.Variant_Classification,'Frame_Shift_Ins'));
M.Start_position=cellfun(@str2num,M.Start_position);
amplicons.gene=[];
amplicons.chr=[];
amplicons.am=[];
for i=1:slength(S)
    
    
k=ismember(M.Hugo_Symbol,S.gene{i});
positions=M.Start_position(k);
positions=sort(positions);
clear am
chr=unique({M.Chromosome{k}});
am(1,1)=positions(1);
    for j=1:length(positions)
    if positions(j)>am(end)+100
        am(end+1,1)=positions(j);
    end
    end
amplicons.gene=vertcat(amplicons.gene,repmat(str2cell(S.gene{i}),size(am)));
amplicons.chr=vertcat(amplicons.chr,repmat(chr,size(am)));
amplicons.am=vertcat(amplicons.am,am);



end


end

function test
% Maf='/Users/amaro/Downloads/EBV_Cluster_FINAL.final_analysis_set.maf';
% sig_genes='/Users/amaro/Downloads/EBV_Cluster_FINAL.sig_genes.txt';
% type='STAD';
% EBVamplicons=MAFtoAmplicons(Maf,sig_genes,type);
% 
% Maf='/Users/amaro/Downloads/MSI_High_Cluster_FINAL.final_analysis_set.maf';
% sig_genes='/Users/amaro/Downloads/MSI_High_Cluster_FINAL.sig_genes.txt';
% MSIamplicons=MAFtoAmplicons(Maf,sig_genes,type)
% 
% Maf='/Users/amaro/Downloads/SCNA_High_Cluster_FINAL.final_analysis_set.maf';
% sig_genes='/Users/amaro/Downloads/SCNA_High_Cluster_FINAL.sig_genes.txt';
% SCNAHigh_amplicons=MAFtoAmplicons(Maf,sig_genes,type)
% 
% Maf='/Users/amaro/Downloads/SCNA_Low_Cluster_FINAL.final_analysis_set.maf';
% sig_genes='/Users/amaro/Downloads/SCNA_Low_Cluster_FINAL.sig_genes.txt';
% SCNALow_amplicons=MAFtoAmplicons(Maf,sig_genes,type)
% 

Maf='~/Projects/STAD/STAD_v8/HyperMaf.af30.ac4.txt';
sig_genes='~/Projects/STAD/STAD_v8/HyperMutsig30af4ac/sig_genes.txt'
type='blah';
Hyperamplicons=MAFtoAmplicons(Maf,sig_genes,type);


Maf='~/Projects/STAD/STAD_v8/NormalNoIntrons.txt';
sig_genes= '~/Projects/STAD/STAD_v8/NormalSigGenesForComut.tsv';
normamplicons=MAFtoAmplicons(Maf,sig_genes,type);

 amplicons=mergeStruct(Hyperamplicons,normamplicons);
   
 
 
 genes=unique(amplicons.gene);
 final_amplicons.gene=[];
 final_amplicons.chr=[];
 final_amplicons.pos=[];
 
 for i=1:length(genes)
    k=find(ismember(amplicons.gene,genes{i})==1);
    positions=sort(amplicons.am(k));
    chr=unique({amplicons.chr{k}});
    clear f_am
    f_am(1,1)=positions(1);
    
    for j=1:length(positions)
        if positions(j)>f_am(end)+100
        f_am(end+1,1)=positions(j);
        end
    end
final_amplicons.gene=vertcat(final_amplicons.gene,repmat(str2cell(genes{i}),size(f_am)));
final_amplicons.chr=vertcat(final_amplicons.chr,repmat(chr,size(f_am)));
final_amplicons.pos=vertcat(final_amplicons.pos,f_am); 
    
 end
 save_struct(final_amplicons,'~/Projects/STAD/STAD_v8/AmpliconsList.txt');
 

end
