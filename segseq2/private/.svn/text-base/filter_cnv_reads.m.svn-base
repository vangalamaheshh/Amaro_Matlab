function [READN,READT]=filter_cnv_reads( IN_READN, IN_READT, cnvFile )

CNV = idx_cnv_reformat(cnvFile);

idxKeepN = ones(length(IN_READN.chr),1);
idxKeepT = ones(length(IN_READT.chr),1);

chrList = 1:max(IN_READN.chr);
fprintf(1,['Filtering CNV in ' cnvFile '\n']);

for c=1:length(chrList)
    chr=chrList(c);
    fprintf(1,'%d..',chr);

    idxChrCNV = find(CNV.chr==chr);
    idxChrN = find(IN_READN.chr==chr);
    idxChrT = find(IN_READT.chr==chr);
    
    posN = IN_READN.pos(idxChrN);
    posT = IN_READT.pos(idxChrT);

    for i = 1:length(idxChrCNV)
	cnvL = CNV.start(idxChrCNV(i));
        cnvR = CNV.end(idxChrCNV(i));
        idxCurrN = find(posN>=cnvL & posN<=cnvR);
        idxCurrT = find(posT>=cnvL & posT<=cnvR);

        idxKeepN(idxChrN(idxCurrN)) = 0;
        idxKeepT(idxChrT(idxCurrT)) = 0;

    end
end

    fprintf(1,'\nNormal reads in CNV:  %.0f\n', length(idxKeepN) - sum(idxKeepN) );
    fprintf(1,'Tumor  reads in CNV:  %.0f\n', length(idxKeepT) - sum(idxKeepT) );

READN.chr = IN_READN.chr(find(idxKeepN==1));
READN.pos = IN_READN.pos(find(idxKeepN==1));
READN.lane = IN_READN.lane(find(idxKeepN==1));

READT.chr = IN_READT.chr(find(idxKeepT==1));
READT.pos = IN_READT.pos(find(idxKeepT==1));
READT.lane = IN_READT.lane(find(idxKeepT==1));
