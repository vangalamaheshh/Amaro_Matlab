function removeCNV( T, germlineCNV );
 
    CNV = idx_cnv_reformat(germlineCNV);

  idxKeepN = ones(length(T.READN.chr),1);
  idxKeepT = ones(length(T.READT.chr),1);

  for ci=1:length(chrs)
    chr=chrs(ci);
    idxChrCNV = find(CNV.chr==chr);
    idxChrN = find(T.READN.chr==chr);
    idxChrT = find(T.READT.chr==chr);
    
    posN = T.READN.pos(idxChrN);
    posT = T.READT.pos(idxChrT);

    for i = 1:length(idxChrCNV)
	cnvL = CNV.start(idxChrCNV(i));
        cnvR = CNV.end(idxChrCNV(i));
        idxCurrN = find(posN>=cnvL & posN<=cnvR);
        idxCurrT = find(posT>=cnvL & posT<=cnvR);

        idxKeepN(idxChrN(idxCurrN)) = 0;
        idxKeepT(idxChrT(idxCurrT)) = 0;
    end
  end
    fprintf(1,'Normal reads in CNV:  %.0f\n', length(idxKeepN) - sum(idxKeepN) );
    fprintf(1,'Tumor  reads in CNV:  %.0f\n', length(idxKeepT) - sum(idxKeepT) );
idxKeepN=logical(idxKeepN);
idxKeepT=logical(idxKeepT);

T.READN.chr = T.READN.chr(idxKeepN));
T.READN.pos = T.READN.pos(idxKeepN);
T.READN.lane = T.READN.lane(idxKeepN);

T.READT.chr = T.READT.chr(idxKeepT);
T.READT.pos = T.READT.pos(idxKeepT);
T.READT.lane = T.READT.lane(idxKeepT);
