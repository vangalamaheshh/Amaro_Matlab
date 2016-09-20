 A.CCF_hat=str2double(A.CCF_hat);
 A.ID=repmat({'NaN'},slength(A),1);
  A.IS_SCNA=str2double(A.IS_SCNA);
  
  A.length=str2double(A.length);
  A.Chromosome=chrom2num(A.Chromosome);
  A.Start_bp=str2double(A.Start_bp);
  A.End_bp=str2double(A.End_bp); 
  A.NA=str2double(A.NA);
  A.NB=str2double(A.NB);
  

lia=( A.CCF_hat>.2 & A.Chromosome==17 & A.Start_bp<7569720 & A.End_bp>7569720 )  &(A.NA==0|A.NB==0);
A.ID(lia,1) = {'del17p'};
% missing GCLL-0001 %some false positives

lia = ( A.CCF_hat>.2 & A.IS_SCNA & A.Chromosome==2 & A.Start_bp>1 & A.End_bp<92719545   & A.length >15000000 ) &(A.NA+A.NB)>2; 
A.ID(lia,1) = {'amp2p'};
% missing CLL-GCLL-0288 (hyper segmented) %some false positives



lia = (A.CCF_hat>.2 & A.Chromosome==6 & A.Start_bp<110014415 & A.End_bp<200998247 & A.length >1000000 )&(A.NA+A.NB)<2;
A.ID(lia,1) = {'del6q21'};


lia = (A.CCF_hat>.2 & A.Chromosome==8 & A.Start_bp<125555423 & A.End_bp>110566132 & A.length >1000000 ) &(A.NA+A.NB)>2;
A.ID(lia,1) = {'amp8q'};


%%del8p

lia = (A.CCF_hat>.2 & A.Chromosome==8 & A.Start_bp>1 & A.End_bp<49887407  & A.length > 1000000 ) & ((A.NA+A.NB)<2|A.NA==0|A.NB==0);
A.ID(lia,1) = {'del8p'};



%%del11q
lia = (A.CCF_hat>.2 & A.Chromosome==11 & A.Start_bp<144898510 & A.End_bp>54974654 & A.length > 1000000  ) &((A.NA+A.NB)==1|A.NA==0|A.NB==0);        
A.ID(lia,1) = {'del11q'};
%CLL-CW106-Tumor-SM-Z14I ?

%tri 12
lia = (A.CCF_hat>.2 & A.Chromosome==12 & A.length>2000000 ) &(A.NA+A.NB)==3;
A.ID(lia,1) = {'tri12'};


%del13q
lia = (A.CCF_hat>.2 & A.Chromosome==13 & A.Start_bp>18887860  ) &((A.NA+A.NB)==1|A.NA==0|A.NB==0);
A.ID(lia,1) = {'del13q'};


%del18p
lia = (A.CCF_hat>.2 & A.Chromosome==18 & A.Start_bp<18835136 & A.length>1083513) &((A.NA+A.NB)==1|A.NA==0|A.NB==0);
A.ID(lia,1) = {'del18p'};

%tri 19
 lia = (A.CCF_hat>.2 & A.Chromosome==19 & A.length>12279373) &(A.NA+A.NB)>2;
 A.ID(lia,1) = {'tri19'}; 
 
 %del20p
 lia = (A.CCF_hat>.2 & A.Chromosome==20 & A.End_bp<46000000 & A.length>2813530) &((A.NA+A.NB)==1|A.NA==0|A.NB==0);
A.ID(lia,1) = {'del20p'};
