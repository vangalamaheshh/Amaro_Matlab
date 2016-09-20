% Amaro Taylor-Weiner Supplementary Code
% #
% #  Allele frequency in African and European population
% #  Tab separated
% #
% #  Gene name - see table for annotation
% #  Location -	      non-coding (----- or UTR),
% #                     synonymous substitution (SYNON)
% #                     non-synonymous substitution (NON-SYNON)
% #  AF - Allele frequency in the African Descent population
% #  EF - Allele frequency in the European Descent population
% #  aa1 -- amino acid (--- if non-coding)
% #  aa2 -- amino acid (--- if non-coding)

SNP_data=load_struct('~/Downloads/PS2_data.dat');
SNP_data.AF=str2double(SNP_data.AF);
SNP_data.EF=str2double(SNP_data.EF);

%Write a program to calculate Tajimas?D and its P-value (described in the notes) Estimate how much variation a population can sustain.
% pi = avg pair wise difference amoung individuals
% pi
% S = number of variable sites
% theta = expectation on pi in neutral population no selection
% theta = S / sum(1./[1:N-1])  N is the number of sequences (people)
% Tajimas D pi (observed) - theta (expected)/ stdev(d)
% caluclate average pairwise difference : http://web.mit.edu/hst.508/HST.508_Biophysics_170/TajimaD_calculations.pdf
% See function TajimasD.m


% part A
NAF=24;NEU=23; % Number of individuals
SNPEU=reorder_struct(SNP_data,SNP_data.EF<1&SNP_data.EF>0);
SNPAF=reorder_struct(SNP_data,SNP_data.AF<1&SNP_data.AF>0);
[D_EU PEU]=TajimasD(NEU,SNPEU,'EF')
[D_AF PAF]=TajimasD(NAF,SNPAF,'AF')



% part B C and D
UTR5=SNP_data.LocationType{12};
UTR3=SNP_data.LocationType{114};
NACLASS=SNP_data.LocationType{1};
SYNON=SNP_data.LocationType{108};
NONSYN=SNP_data.LocationType{110};
TRUNC=SNP_data.LocationType{4654};
FS=SNP_data.LocationType{4655};
IFS=SNP_data.LocationType{25128};
SPLICE=SNP_data.LocationType{26924};


% European Data Set
CodingMuts=reorder_struct(SNP_data,~ismember(SNP_data.LocationType,{UTR3,UTR5,NACLASS})&SNP_data.EF<1&SNP_data.EF>0);
NonCodingMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{UTR3,UTR5,NACLASS})&SNP_data.EF<1&SNP_data.EF>0);
[D_CodeEU PCodeEU]=TajimasD(NEU,CodingMuts,'EF')
[D_NonCodeEU P_NoncodeEU]=TajimasD(NEU,NonCodingMuts,'EF')

SynMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{SYNON})&SNP_data.EF<1&SNP_data.EF>0);
NonSynMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{TRUNC,NONSYN,FS,IFS,SPLICE})&SNP_data.EF<1&SNP_data.EF>0);
[D_SynEU, PSynEU]=TajimasD(NEU,SynMuts,'EF')
[D_NonSynEU, PNonSynEU]=TajimasD(NEU,NonSynMuts,'EF')

CodingMuts=reorder_struct(SNP_data,~ismember(SNP_data.LocationType,{UTR3,UTR5,NACLASS})&SNP_data.AF<1&SNP_data.AF>0);
NonCodingMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{UTR3,UTR5,NACLASS})&SNP_data.AF<1&SNP_data.AF>0);
[D_CodeAF PCodeAF]=TajimasD(NAF,CodingMuts,'AF')
[D_NonCodeAF PNonCodeAF]=TajimasD(NAF,NonCodingMuts,'AF')

SynMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{SYNON})&SNP_data.AF<1&SNP_data.AF>0);
NonSynMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{TRUNC,NONSYN,FS,IFS,SPLICE})&SNP_data.AF<1&SNP_data.AF>0);
[D_SynAF PSynAF]=TajimasD(NAF,SynMuts,'AF')
[D_NonSynAF PNonSynAF]=TajimasD(NAF,NonSynMuts,'AF')

CodingSynMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{SYNON})&SNP_data.EF<1&SNP_data.EF>0);
NonCodingMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{UTR3,UTR5,NACLASS})&SNP_data.EF<1&SNP_data.EF>0);


D_SynEU=TajimasD(NEU,CodingSynMuts,'EF')
D_NonCodingEU=TajimasD(NEU,NonCodingMuts,'EF')

UTRMUTS=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{UTR3,UTR5})&SNP_data.EF<1&SNP_data.EF>0);
SynMuts=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{SYNON})&SNP_data.EF<1&SNP_data.EF>0);

[DUTREU PEUUTR]=TajimasD(NEU,UTRMUTS,'EF')

UTRMUTS=reorder_struct(SNP_data,ismember(SNP_data.LocationType,{UTR3,UTR5})&SNP_data.AF<1&SNP_data.AF>0);
[DUTREU PEUUTR]=TajimasD(NAF,UTRMUTS,'AF')



order={'ALA','ARG', 'ASN','ASP','CYS','GLN','GLU','GLY','HIS' ,'ILE' ,'LEU' ,'LYS','MET','PHE' ,'PRO', 'SER','THR' ,'TRP','TYR','VAL'};
BLOSUM_MATRIX=(blosum(50));

for i=1:slength(SNP_data_nononcode)
    
    x=find(ismember(order,SNP_data_nononcode.aa1{i}));
    y=find(ismember(order,SNP_data_nononcode.aa2{i}));
    SNP_data_nononcode.blosum(i,1)=BLOSUM_MATRIX(x,y);
    
end

Hydrophobic={'ALA'
'ILE'
'LEU'
'PHE'
'VAL'
'PRO' 
'GLY'
};
for i=1:slength(SNP_data_nononcode)
    if (ismember(SNP_data_nononcode.aa1{i},Hydrophobic)&&ismember(SNP_data_nononcode.aa2{i},Hydrophobic))||...
            (~ismember(SNP_data_nononcode.aa1{i},Hydrophobic)&&~ismember(SNP_data_nononcode.aa2{i},Hydrophobic))
        SNP_data_nononcode.hydro(i,1)=0;
    else
        
    SNP_data_nononcode.hydro(i,1)=1;
    end
    
end

