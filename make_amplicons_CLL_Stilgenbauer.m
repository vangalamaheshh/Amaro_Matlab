% always add NOTCH1 Indel
amplicons.Start(1)=139390649;
amplicons.End(1)=139390649;
amplicons.Chr{1}='9';
amplicons.Gene{1}='NOTCH1';
amplicons.Sample{1}='EverySample';
% GenesToValidateOnlyIfSalvage={'SF3B1'
% 'ATM'
% 'TP53'
% 'POT1'
% 'NOTCH1'
% 'XPO1'
% 'BIRC3'
% 'BRAF'
% 'EGR2'
% 'MYD88'
% 'KRAS'
% 'DDX3X'
% 'NRAS'
% 'SAMHD1'
% 'FBXW7'
% 'ZMYM3'
% 'PTPN11'
% 'MED12'
% 'BCOR'
% };

CodingClass=get_coding_class_muts();

%Redundant_Samples={'CLL-GCLL-0233-Tumor-SM-4DP99';'CLL-GCLL-0049-Tumor-SM-41JND';'CLL-GCLL-0157-Tumor-SM-41JYO';'CLL-GCLL-0164-Tumor-SM-41JYV';'CLL-GCLL-0022-Tumor-SM-41JML'};
%285,233,166,49,157,164,68,22  are on both chips keep 68 and 166

M=load_struct('/Users/amaro/Documents/StilgenbauerLongChips/chip3.txt');
%MStil=load_struct('/Volumes/xchip_cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.maf');
%M=load_struct('/xchip/cga_home/amaro/CLL/StilgenbauerMafFreeze2.0/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.maf');

   %DFCI_Samples={'CLL-CR-Tumor';'CLL-CW106-Tumor-SM-Z14I';'CLL-CW107-Tumor-SM-Z133';'CLL-CW108-Tumor-SM-Z135';'CLL-CW109-Tumor-SM-Z137';'CLL-CW11-Tumor-SM-UFNK';'CLL-CW113-Tumor-SM-19DUP';'CLL-CW117-Tumor-SM-1EVH7';'CLL-CW12-Tumor-SM-UFNM';'CLL-CW127-Tumor-SM-19DV2';'CLL-CW130-Tumor-SM-1EVHC';'CLL-CW159-Tumor-SM-19DW4';'CLL-CW161-Tumor-SM-1EVI1';'CLL-CW163-Tumor-SM-19DW8';'CLL-CW165-Tumor-SM-19DWA';'CLL-CW207-Tumor-SM-1EVRK';'CLL-CW210-Tumor-SM-1EVRN';'CLL-CW213-Tumor-SM-1EVRQ';'CLL-CW216-Tumor-SM-1EVRT';'CLL-CW218-Tumor-SM-1EVRV';'CLL-CW224-Tumor-SM-1EVS2';'CLL-CW227-Tumor-SM-1EVS5';'CLL-CW23-Tumor-SM-UFO9';'CLL-CW230-Tumor-SM-1EVS8';'CLL-CW236-Tumor-SM-1EVSE';'CLL-CW24-Tumor-SM-UFOB';'CLL-CW244-Tumor-SM-1EVSM';'CLL-CW247-Tumor-SM-1EVSP';'CLL-CW248-Tumor-SM-1EVSQ';'CLL-CW27-Tumor-SM-UFOH';'CLL-CW32-Tumor-SM-UFOR';'CLL-CW45-Tumor-SM-24ABJ';'CLL-CW49-Tumor-SM-UFPQ';'CLL-CW50-Tumor-SM-UFPS';'CLL-CW58-Tumor-SM-UFRC';'CLL-CW64-Tumor-SM-UFRO';'CLL-CW73-Tumor-SM-UFS7';'CLL-CW82-Tumor-SM-UFSP';'CLL-CW89-Tumor-SM-UFT4';'CLL-CW_194-Tumor-SM-19DX4';'CLL-CW_202-Tumor-SM-19DXC';'CLL-DK-TN-Tumor';'CW149-TN-Tumor';'CW85-PRE-Tumor'};
   % M=reorder_struct(M,ismember(M.Tumor_Sample_Barcode,DFCI_Samples));
%M=reorder_struct(M,ismember(M.Tumor_Sample_Barcode,MStil.Tumor_Sample_Barcode));
sig_genes=load_struct('/Volumes/xchip_cga_home/amaro//CLL/StilgenbauerMafFreeze2.0/ICGC_StilgenbauerForComut1.16/ICGC_Stilgenbauer.aggregate.NOIGV_WITHTTKandNOTCH1_PoNCut4_WithSalvage.sig_genes_.txt');
M.Start_position=str2double(M.Start_position);
sig_genes.q=str2double(sig_genes.q);
sig_genes=reorder_struct(sig_genes,sig_genes.q<=.1);
sig_genes=reorder_struct(sig_genes,~ismember(sig_genes.gene,{'ADAM30';'NHS';'TTK'}));
M=reorder_struct(M,ismember(M.Hugo_Symbol,sig_genes.gene));
M=reorder_struct(M,~ismember(M.Hugo_Symbol,GenesToValidateOnlyIfSalvage)|ismember(M.PoN_Status,'SALVAGE'));
%M=reorder_struct(M,ismember(M.Variant_Classification,CodingClass));
%M=reorder_struct(M,~ismember(M.Tumor_Sample_Barcode,Redundant_Samples));
l_amp=1;

for i=1:slength(M)
    lc=find(ismember(amplicons.Chr,M.Chromosome{i}));
    flipc=0;
    if isempty(lc)
        amplicons.Start(l_amp+1)=M.Start_position(i)-75;
        amplicons.End(l_amp+1)=M.Start_position(i)+75;
        amplicons.Chr{l_amp+1}=M.Chromosome{i};
        amplicons.Gene{l_amp+1}=M.Hugo_Symbol{i};
        amplicons.Sample{l_amp+1}=M.Tumor_Sample_Barcode{i};
        l_amp=l_amp+1;
    else
        for j=1:length(lc)
        if amplicons.Start(lc(j))<=M.Start_position(i)&&M.Start_position(i)<=amplicons.End(lc(j))
            amplicons.Sample{lc(j)}=strcat(amplicons.Sample{lc(j)},',',M.Tumor_Sample_Barcode{i});
            flipc=1;
        end
        end
        if flipc==0
        amplicons.Start(l_amp+1)=M.Start_position(i)-75;
        amplicons.End(l_amp+1)=M.Start_position(i)+75;
        amplicons.Chr{l_amp+1}=M.Chromosome{i};
        amplicons.Gene{l_amp+1}=M.Hugo_Symbol{i};
        amplicons.Sample{l_amp+1}=M.Tumor_Sample_Barcode{i};
        l_amp=l_amp+1;
        end
    end

end

  
amplicons.Start=amplicons.Start';
amplicons.End=amplicons.End';
amplicons.Chr=amplicons.Chr';
amplicons.Gene=amplicons.Gene';
amplicons.Sample=amplicons.Sample';

