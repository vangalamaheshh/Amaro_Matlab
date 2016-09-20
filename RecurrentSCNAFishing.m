%SC=load_struct('/Volumes/xchip_cga_home/amaro/CLL/ABSOLUTE/AggregateCCF_hat.Stil.DFCI.ICGC.seg');
 Chroms = {'2', '6', '8', '11', '12', '13', '17', '18', '19', '20'};
 Lia = ismember(A.Chromosome, Chroms);
 A = reorder_struct(A, Lia);
 A.CCF_hat=str2double(A.CCF_hat);
 A.ID=repmat({'NaN'},slength(A),1);
  A.IS_SCNA=str2double(A.IS_SCNA);
  
  A.length=str2double(A.length);
  A.Chromosome=chrom2num(A.Chromosome);
  A.Start_bp=str2double(A.Start_bp);
  A.End_bp=str2double(A.End_bp); 
  A.NA=str2double(A.NA);
  A.NB=str2double(A.NB);
  
  
 loc = 61767418; 
 Lia = (A.IS_SCNA & A.Chromosome==2 & A.Start_bp<loc & A.End_bp>loc & A.length >1000000 ) &(A.NA+A.NB)==3;
 blacklist = {
    'CLL-CW139-TP-NT-SM-19DVB-SM-1EVUY'
    'CLL-CW163-TP-NT-SM-19DW8-SM-19E7L'
    'CLL-CW53-TP-NT-SM-UFPY-SM-UFPZ'
    'CLL-DK-TN-TB-NBC--'
    'CLL-GCLL-0240-TP-NT-SM-4DP9G-SM-4DOJR'
    'CLL-GCLL-0138-TP-NT-SM-41JY5-SM-41JSN'
        };
idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
 A.ID(Lia,1) = {'amp2p'};

 
%%del6q21
loc = 110014415;
Lia = (A.IS_SCNA & A.Chromosome==6 & A.Start_bp<loc & A.End_bp>loc & A.length >1000000 )&(A.NA+A.NB)==1;
blacklist = {
    'ICGC_CLL-166-TB-NBC--'
    'ICGC_CLL-642-TB-NBC--'
    'CLL-CW124-TP-NT-SM-24ABT-SM-24ABU'
    'CLL-GCLL-0246-TP-NT-SM-4DP9M-SM-4DOJX'
    'CLL-GCLL-0158-TP-NT-SM-41JYP-SM-41JT8'
    'CLL-CW58-TP-NT-SM-UFRC-SM-UFRD'
        };
        

idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'del6q21'};

 
%%amp8q
loc = 128746315;
Lia = (A.IS_SCNA & A.Chromosome==8 & A.Start_bp<loc & A.End_bp>loc & A.length >1000000 ) &(A.NA==2 | A.NB==2);
blacklist = {
    'CLL-GCLL-0127-TP-NT-SM-41JXT-SM-41JSC'
    'CLL-GCLL-0246-TP-NT-SM-4DP9M-SM-4DOJX'
    'CLL-CR-TB-NBC--'
    'CLL-CW41-TP-NT-SM-UFPA-SM-UFPB'
    'CLL-CW48-TP-NT-SM-UFPO-SM-UFPP'
    'CLL-DK-TN-TB-NBC--'
       };
        
        

idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'amp8q'};

%%del8p
zone_end = 43124113;
zone_start = 1;
Lia = (A.IS_SCNA & A.Chromosome==8 & A.Start_bp>zone_start & A.End_bp<zone_end & A.n_probes>50  ) & (A.minNA+A.minNB==1 | A.minNA+A.minNB==0);
blacklist = {
    'CLL-GCLL-0127-TP-NT-SM-41JXT-SM-41JSC'
    'CLL-GCLL-0246-TP-NT-SM-4DP9M-SM-4DOJX'
    'CLL-CR-TB-NBC--'
    'CLL-CW41-TP-NT-SM-UFPA-SM-UFPB'
    'CLL-CW48-TP-NT-SM-UFPO-SM-UFPP'
    'CLL-DK-TN-TB-NBC--'
       };
        
        

idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'del8p'};




%%del11q
loc = 108091559;
Lia = (A.IS_SCNA & A.Chromosome==11 & A.Start_bp<loc & A.End_bp>loc  ) &(A.NA+A.NB)==1;
blacklist = {
    'CLL-GCLL-0127-TP-NT-SM-41JXT-SM-41JSC'
    'CLL-GCLL-0246-TP-NT-SM-4DP9M-SM-4DOJX'
    'CLL-CR-TB-NBC--'
    'CLL-CW41-TP-NT-SM-UFPA-SM-UFPB'
    'CLL-CW48-TP-NT-SM-UFPO-SM-UFPP'
    'CLL-DK-TN-TB-NBC--'
       };
        
        

idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'del11q'};


%%tri12
loc = 108091559;

Lia = (A.IS_SCNA & A.Chromosome==12 & A.Start_bp<loc & A.End_bp>loc & A.length>2000000 ) &(A.NA+A.NB)==3;

blacklist = {
'CLL-GCLL-0138-TP-NT-SM-41JY5-SM-41JSN'
'CLL-CR-TB-NBC--'
'CLL-CW48-TP-NT-SM-UFPO-SM-UFPP'
'CLL-DK-TN-TB-NBC--'
'CLL-JP-TB-NBC--'
        };
        
idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'tri12'};



%%del13q
loc = 50400000;
Lia = (A.IS_SCNA & A.Chromosome==13 & A.Start_bp<loc & A.End_bp>loc)  &((A.NA+A.NB)==1);  

blacklist = {
'CLL-GCLL-0178-TP-NT-SM-41JZA-SM-41JTS'
        };

    idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'del13q'};

%del17p

%zone A.Start_bp>=1&A.End_bp<21051758
loc = 7592868;
Lia = (A.IS_SCNA & A.Chromosome==17 & A.Start_bp<loc & A.End_bp>loc & A.length>300000)|( A.CCF_hat>.2 & A.Chromosome==17 & A.Start_bp<7569720 & A.End_bp>7569720 )  &(A.NA==0|A.NB==0);
lia=( A.CCF_hat>.06 & A.Chromosome==17 & A.Start_bp<7569720 & A.End_bp>7569720 )  &(A.NA==0|A.NB==0);

blacklist = {
'CLL-GCLL-0012-TP-NT-SM-41JMB-SM-41QBS'
'CLL-GCLL-0020-TP-NT-SM-41JMJ-SM-41QC1'
'CLL-GCLL-0083-TP-NT-SM-41JOC-SM-41QDT'
'CLL-GCLL-0107-TP-NT-SM-41JX9-SM-41JRR'
'CLL-CW251-TP-NT-SM-1EVST-SM-1EVU7'
 'CLL-CW163-TP-NT-SM-19DW8-SM-19E7L'
 'CLL-CW58-TP-NT-SM-UFRC-SM-UFRD'
        };
        
   
    
idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'del17p'};

whitelist = {
'CLL-GCLL-0239-TP-NT-SM-4DP9F-SM-4DOJQ'
'CLL-GCLL-0001-TP-NT-SM-41JLZ-SM-41QBH'
        };
        
idx = ismember(A.sample, whitelist) & &A.Chromosome==17&A.Start_bp<7569720 & A.End_bp>7569720;
Lia(idx) = 1;
A.ID(Lia,1) = {'del17p'};





%%%del18p
loc = 4056835;
Lia = (A.IS_SCNA & A.Chromosome==18 & A.Start_bp<loc & A.End_bp>loc & A.length>1000000 ) &(A.NA+A.NB)==1;

blacklist = {
'CLL-GCLL-0178-TP-NT-SM-41JZA-SM-41JTS'
'CLL-CW166-TP-NT-SM-19DWB-SM-19E7O'
        };
idx = ismember(A.sample, blacklist);
Lia(idx) = 0;

A.ID(Lia,1) = {'del18p'};


%tri 19
loc = 10056835;
Lia = (A.IS_SCNA & A.Chromosome==19 & A.Start_bp<loc & A.End_bp>loc  & A.length>1000000 ) |(A.IS_SCNA & A.Chromosome==19 & A.Start_bp<20056835 & A.End_bp>20056835  & A.length>1000000) &(A.NA==2 | A.NB==2 );
blacklist = {
 'CLL-GCLL-0045-TP-NT-SM-41JN9-SM-41QCQ'
    'CLL-GCLL-0127-TP-NT-SM-41JXT-SM-41JSC'
    'CLL-GCLL-0138-TP-NT-SM-41JY5-SM-41JSN'
    'CLL-GCLL-0037-TP-NT-SM-41JN1-SM-41QCI'
    'CLL-GCLL-0150-TP-NT-SM-41JYH-SM-41JSZ'
     'CLL-CW230-TP-NT-SM-1EVS8-SM-1EVTL'
    'CLL-CW83-TP-NT-SM-UFSR-SM-UFSS'
        'CLL-CW103-TP-NT-SM-UFTW-SM-UFTX'
    'CLL-CW105-TP-NWGA-SM-Z12Y-SM-Z12Z'
    'CLL-CW109-TP-NT-SM-Z137-SM-Z138'
    'CLL-CW138-TP-NT-SM-19DVA-SM-19E6U'
    'CLL-CW50-TP-NT-SM-UFPS-SM-UFPT'
    'CLL-CW53-TP-NT-SM-UFPY-SM-UFPZ'
    'CLL-CW_176-TP-NT-SM-19DWL-SM-1EVVU'
    'ICGC_CLL-282-TB-NBC--'

           };


idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'tri19'};

%%del20p


loc = 20056835;
Lia = (A.IS_SCNA & A.Chromosome==20 & A.Start_bp<loc & A.End_bp>loc  & A.length>1000000) &(A.NA+A.NB)==1;

blacklist = {
  'CLL-CW58-TP-NT-SM-UFRC-SM-UFRD'
    'CLL-CW_179-TP-NT-SM-19DWO-SM-19E82'
    'CLL-GCLL-0001-TP-NT-SM-41JLZ-SM-41QBH'
    'CLL-GCLL-0012-TP-NT-SM-41JMB-SM-41QBS'
    'CLL-GCLL-0020-TP-NT-SM-41JMJ-SM-41QC1'
    'CLL-GCLL-0083-TP-NT-SM-41JOC-SM-41QDT'
    'CLL-GCLL-0107-TP-NT-SM-41JX9-SM-41JRR'
           };
       
idx = ismember(A.sample, blacklist);
Lia(idx) = 0;
A.ID(Lia,1) = {'del20p'};
